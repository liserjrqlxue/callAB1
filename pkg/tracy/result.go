package tracy

import (
	"fmt"
	"io"
	"log/slog"
	"strings"

	"github.com/liserjrqlxue/goUtil/fmtUtil"
)

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

	SV bool `json:"sv"`
}

func (v *Variant) String() string {
	return fmt.Sprintf("%s\t%d\t%s\t%s\t%s\t%d\t%s\t%s\t%s\t%d\t%d\t%d\t%d",
		v.Chr, v.Pos, v.ID, v.Ref, v.Alt, v.Qual, v.Filter, v.Type, v.Genotype, v.Basepos, v.Signalpos, v.Xrange[0], v.Xrange[1])
}

type Variants struct {
	Columns []string `json:"columns"`
	Rows    [][]any  `json:"rows"`
	Xranges [][2]int `json:"xranges"`

	Variants []*Variant `json:"variants"`

	HetCount int
	HomCount int
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

		if varint.Genotype == "het." {
			v.HetCount++
		} else {
			v.HomCount++
		}

		if varint.Type == "Complex" || (varint.Type == "SNV" && len(varint.Alt) >= SVThreshold) || len(varint.Alt)-len(varint.Ref) >= SVThreshold || len(varint.Ref)-len(varint.Alt) >= SVThreshold {
			varint.SV = true
		}
	}
}

func (v Variants) String() string {
	var lines []string
	for _, v := range v.Variants {
		lines = append(lines, v.String())
	}
	return strings.Join(lines, "\n")
}

type Result struct {
	Meta *Meta `json:"meta"`

	Pos   []int //`json:"pos"`
	PeakA []int //`json:"peakA"`
	PeakC []int //`json:"peakC"`
	PeakG []int //`json:"peakG"`
	PeakT []int //`json:"peakT"`

	BasecallPos  []int `json:"basecallPos"`
	BasecallQual []int `json:"basecallQual"`

	Basecalls map[string]string `json:"basecalls"`

	PrimarySeq   string `json:"primarySeq"`
	SecondarySeq string `json:"secondarySeq"`

	ChartConfig map[string]any `json:"chartConfig"`

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

	AlignResult *AlignResult `json:"alignResult"`

	Status string `json:"status"`
	Pass   bool   `json:"pass"`
}

type Trace struct {
	TraceFileName string `json:"traceFileName"`

	LeadingGaps  int `json:"leadingGaps"`
	TrailingGaps int `json:"trailingGaps"`

	PeakA []int `json:"peakA"`
	PeakC []int `json:"peakC"`
	PeakG []int `json:"peakG"`
	PeakT []int `json:"peakT"`

	BasecallPos  []int `json:"basecallPos"`
	BasecallQual []int `json:"basecallQual"`

	Basecalls map[string]string `json:"basecalls"`
}

type AlignResult struct {
	GappedTrace *Trace `json:"gappedTrace"`

	Refchr   string `json:"refchr"`
	Refpos   int    `json:"refpos"`
	Altalign string `json:"altalign"`
	Refalign string `json:"refalign"`
	Forward  int    `json:"forward"`

	BoundMatch      int     `json:"boundMatch"`
	BoundMatchRatio float64 `json:"boundMatchRatio"`
	BoundStatus     []int   `json:"boundStatus"`

	MatchRegion [2]int `json:match_region`
}

func (ar *AlignResult) CalAlign() {
	if len(ar.Refalign) != len(ar.Altalign) {
		panic("len(ar.Refalign)!=len(ar.Altalign)")
	}
	if ar.Altalign[0] == '-' || ar.Altalign[len(ar.Altalign)-1] == '=' {
		panic("ar.Altalign start/stop with -")
	}
	var (
		match     = 0
		refLength = 0
		altLength = 0
		// altInsert = 0
		// altDel    = 0
		// 0:match 1:insert 2:del 3:mismatch
		currentStatus = 0
		status        []int

		leftCut  = 0
		rightCut = 0
	)
	alignLength := len(ar.Refalign)
	for i := range ar.Refalign {
		refN := ar.Refalign[i]
		if refN == '-' {
			leftCut++
		} else {
			break
		}
	}

	for i := range ar.Refalign {
		refN := ar.Refalign[alignLength-1-i]
		if refN == '-' {
			rightCut++
		} else {
			break
		}
	}

	for i := range ar.Refalign {
		refN := ar.Refalign[i]
		altN := ar.Altalign[i]

		if refN == altN {
			if refN == '-' {
				panic("refN=='-'&&altN=='-'")
			}
			currentStatus = 0
			refLength++
			altLength++
			match++
		} else if refN == '-' {
			currentStatus = 1
			altLength++
		} else if altN == '-' {
			currentStatus = 2
			refLength++
		} else {
			currentStatus = 3
			refLength++
			altLength++
		}
		status = append(status, currentStatus)
	}

	// match region 0-base
	// forward:
	//   refpos                  refLength + refpos-1
	//   leftCut+1               altLength-rightCut
	// reverse:
	//   refLenth + refpos-1     refpos
	//   leftCut+1               altLength-rightCut
	// fmt.Printf("Match Region: forward:%d:Ref[%d,%d] Alt[%d,%d]\n", ar.Forward, ar.Refpos, ar.Refpos+refLength-1, leftCut+1, altLength-rightCut)
	leftOffset := max(0, Trim-leftCut)
	rightOffset := max(0, altLength-rightCut-MaxLength)
	if ar.Forward == 1 {
		// fmt.Printf("Match Region: forward:%d:Ref[%d,%d] Alt[%d,%d]\n", ar.Forward, ar.Refpos+leftOffset, ar.Refpos+refLength-1-rightOffset, leftCut+1+leftOffset, altLength-rightCut-rightOffset)
		ar.MatchRegion = [2]int{ar.Refpos + leftOffset, ar.Refpos + refLength - 1 - rightOffset}
	} else {
		// fmt.Printf("Match Region: forward:%d:Ref[%d,%d] Alt[%d,%d]\n", ar.Forward, ar.Refpos+rightOffset, ar.Refpos+refLength-1-leftOffset, leftCut+1+leftOffset, altLength-rightCut-rightOffset)
		ar.MatchRegion = [2]int{ar.Refpos + rightOffset, ar.Refpos + refLength - 1 - leftOffset}
	}
	slog.Debug("CalAlign", "match", match, "refLength", refLength, "altLength", altLength, "altRatio", float64(match)/float64(altLength), "status", status)

	var altPos = 0
	var refPos = 0
	var altBound = [2]int{Trim, min(altLength-Trim, MaxLength)}
	var refBound = [2]int{40, refLength - 40}
	var boundStatus []int
	var newBound [4]int
	var boundMatch = 0
	for _, s := range status {
		if s != 1 {
			refPos++
		}
		if s != 2 {
			altPos++
		}
		if altPos >= altBound[0] && altPos <= altBound[1] && refPos >= refBound[0] && refPos <= refBound[1] {
			boundStatus = append(boundStatus, s)
			if len(boundStatus) == 1 {
				newBound[0] = refPos
				newBound[2] = altPos
			}
			newBound[1] = refPos
			newBound[3] = altPos
			if s == 0 {
				boundMatch++
			}
		}
	}
	ar.BoundMatch = boundMatch
	ar.BoundMatchRatio = float64(boundMatch) / float64(len(boundStatus))
	ar.BoundStatus = boundStatus

	slog.Debug("Bound", "refBound", refBound, "altBound", altBound, "newBound", newBound, "boundMatch", boundMatch, "boundMatchRatio", ar.BoundMatchRatio, "boundStatus", ar.BoundStatus)
}

func (ar *AlignResult) Summary(summary io.Writer) {
	fmtUtil.Fprintf(summary, "BoundMatch: %d, BoundLength: %d, BoundMatchRatio: %f\n", ar.BoundMatch, len(ar.BoundStatus), ar.BoundMatchRatio)
}
