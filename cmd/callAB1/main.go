package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"log"
	"log/slog"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"

	"github.com/liserjrqlxue/goUtil/fmtUtil"
	"github.com/liserjrqlxue/goUtil/osUtil"
	"github.com/liserjrqlxue/goUtil/simpleUtil"
)

// flag
var (
	tracy = flag.String(
		"tracy",
		"tracy",
		"https://github.com/gear-genomics/tracy binary path",
	)
	input = flag.String(
		"input",
		"input",
		"input ab1 path",
	)
	ref = flag.String(
		"ref",
		"",
		"ref fasta path",
	)
	prefix = flag.String(
		"prefix",
		"",
		"output prefix",
	)
)

func main() {
	flag.Parse()
	if *input == "" {
		flag.PrintDefaults()
		log.Fatal("input is required")
	}
	if *ref == "" {
		flag.PrintDefaults()
		log.Fatal("ref is required")
	}
	if *prefix == "" {
		*prefix = strings.TrimSuffix(filepath.Base(*input), filepath.Ext(*input))
	}

	var left, right = 0, 0
	var resultJson []byte
	var tracyResult *tracyResult

	simpleUtil.CheckErr(basecall(*tracy, *input, *prefix, left, right))
	resultJson = simpleUtil.HandleError(os.ReadFile(*prefix + ".basecall.json"))

	if err := json.Unmarshal(resultJson, &tracyResult); err != nil {
		log.Fatalf("json.Unmarshal failed: %v\n\tResult:%s", err, resultJson)
	}

	// trim
	left = 60
	right = max(60, len(tracyResult.BasecallPos)-700)
	simpleUtil.CheckErr(decompose(*tracy, *ref, *input, *prefix, left, right))
	resultJson = simpleUtil.HandleError(os.ReadFile(*prefix + ".decompose.json"))

	if err := json.Unmarshal(resultJson, &tracyResult); err != nil {
		log.Fatalf("json.Unmarshal failed: %v\n\tResult:%s", err, resultJson)
	}

	tracyResult.Variants.CalVariants()
	out := osUtil.Create(*prefix + ".decompose.variants.txt")
	defer simpleUtil.DeferClose(out)
	// output title
	fmtUtil.FprintStringArray(out, append(tracyResult.Variants.Columns, "xrange1", "xrange2"), "\t")
	fmtUtil.Fprintln(out, tracyResult.Variants.String())

}

func basecall(tracy, input, prefix string, left, right int) error {
	var args = []string{
		"basecall",
		"-o", prefix + ".basecall.json",
		"-q", strconv.Itoa(left),
		"-u", strconv.Itoa(right),
		input,
	}
	cmd := exec.Command(tracy, args...)
	slog.Info("Basecall", "CMD", cmd)
	cmd.Stderr = os.Stderr
	cmd.Stdout = os.Stdout
	return cmd.Run()
}

func decompose(tracy, ref, input, prefix string, left, right int) error {
	var args = []string{
		"decompose",
		"-r", ref,
		"-o", prefix + ".decompose",
		"-q", strconv.Itoa(left),
		"-u", strconv.Itoa(right),
		"-v",
		input,
	}
	cmd := exec.Command(tracy, args...)
	slog.Info("Decompose", "CMD", cmd)
	cmd.Stderr = os.Stderr
	cmd.Stdout = os.Stdout
	return cmd.Run()
}

// "meta": {"program": "tracy", "version": "0.7.8", "arguments": {"trimLeft": 50, "trimRight": 50, "pratio": 0.33, "genome": "BGA_14A2.fa", "input": "BGA-14A2-TY1-8-T7.ab1"}},
type Meta struct {
	Program   string                 `json:"program"`
	Version   string                 `json:"version"`
	Arguments map[string]interface{} `json:"arguments"`
}

type Variant struct {
	Chr       string `json:"chr"`
	Pos       int    `json:"pos"`
	ID        string `json:"id"`
	Ref       string `json:"ref"`
	Alt       string `json:"alt"`
	Qual      int    `json:"qual"`
	Filter    string `json:"filter"`
	Type      string `json:"type"`
	Genotype  string `json:"genotype"`
	Basepos   int    `json:"basepos"`
	Signalpos int    `json:"signalpos"`

	Xrange [2]int `json:"xrange"`
}

func (v *Variant) String() string {
	return fmt.Sprintf("%s\t%d\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d",
		v.Chr, v.Pos, v.ID, v.Ref, v.Qual, v.Filter, v.Type, v.Genotype, v.Alt, v.Basepos, v.Signalpos, v.Xrange[0], v.Xrange[1])
}

type Variants struct {
	Columns []string        `json:"columns"`
	Rows    [][]interface{} `json:"rows"`
	Xranges [][2]int        `json:"xranges"`

	Variants []*Variant `json:"variants"`
}

func (v *Variants) CalVariants() {
	for i := range v.Rows {
		row := v.Rows[i]
		xrange := v.Xranges[i]
		data := make(map[string]interface{})
		for j := range v.Columns {
			data[v.Columns[j]] = row[j]
		}

		varint := &Variant{
			Chr:       data["chr"].(string),
			Pos:       int(data["pos"].(float64)),
			ID:        data["id"].(string),
			Ref:       data["ref"].(string),
			Alt:       data["alt"].(string),
			Qual:      int(data["qual"].(float64)),
			Filter:    data["filter"].(string),
			Type:      data["type"].(string),
			Genotype:  data["genotype"].(string),
			Basepos:   int(data["basepos"].(float64)),
			Signalpos: int(data["signalpos"].(float64)),

			Xrange: xrange,
		}
		v.Variants = append(v.Variants, varint)
	}
}

func (v Variants) String() string {
	var lines []string
	for _, v := range v.Variants {
		lines = append(lines, v.String())
	}
	return strings.Join(lines, "\n")
}

type tracyResult struct {
	Meta *Meta `json:"meta"`

	Pos   []int `json:"pos"`
	PeakA []int `json:"peakA"`
	PeakC []int `json:"peakC"`
	PeakG []int `json:"peakG"`
	PeakT []int `json:"peakT"`

	BasecallPos  []int             `json:"basecallPos"`
	BasecallQual []int             `json:"basecallQual"`
	Basecalls    map[string]string `json:"basecalls"`
	PrimarySeq   string            `json:"primarySeq"`
	SecondarySeq string            `json:"secondarySeq"`

	ChartConfig map[string]interface{} `json:"chartConfig"`

	Ref1chr     string `json:"ref1chr"`
	Ref1pos     int    `json:"ref1pos"`
	Alt1align   string `json:"alt1align"`
	Ref1align   string `json:"ref1align"`
	Ref1forward int    `json:"ref1forward"`
	Align1score int    `json:"align1score"`

	Ref2chr     string `json:"ref2chr"`
	Ref2pos     int    `json:"ref2pos"`
	Alt2align   string `json:"alt2align"`
	Ref2align   string `json:"ref2align"`
	Ref2forward int    `json:"ref2forward"`
	Align2score int    `json:"align2score"`

	Allele1fraction float64 `json:"allele1fraction"`
	Allele1align    string  `json:"allele1align"`

	Allele2fraction float64 `json:"allele2fraction"`
	Allele2align    string  `json:"allele2align"`

	Align3score int `json:"align3score"`

	Hetindel int `json:"hetindel"`

	Decompositon map[string][]int `json:"decompositon"`

	Variants *Variants `json:"variants"`
}
