package main

import (
	"flag"
	"log"
	"log/slog"
	"os"
	"path/filepath"

	"callAB1/pkg/tracy"

	"github.com/liserjrqlxue/goUtil/fmtUtil"
	"github.com/liserjrqlxue/goUtil/osUtil"
	"github.com/liserjrqlxue/goUtil/simpleUtil"
	"github.com/liserjrqlxue/goUtil/stringsUtil"
	"github.com/xuri/excelize/v2"
)

// flag
var (
	input = flag.String(
		"i",
		"",
		"input 自合.xlsx",
	)
	outputDir = flag.String(
		"o",
		"",
		"output dir",
	)
	sangerDir = flag.String(
		"s",
		"",
		"dir of sanger",
	)
	bin = flag.String(
		"tracy",
		"tracy",
		"https://github.com/gear-genomics/tracy binary path",
	)
)

type Seq struct {
	ID    string `json:"id"`
	Seq   string `json:"seq"`
	Start int    `json:"start"`
	End   int    `json:"end"`

	RefID    string `json:"refID"`
	RefSeq   *Seq   `json:"refSeq"`
	RefStart int    `json:"refStart"`
	RefEnd   int    `json:"refEnd"`

	SubSeq []*Seq `json:"subSeq"`
}

func (s *Seq) CreateSub(id string, start, end int) *Seq {
	var subSeq = &Seq{
		ID:       id,
		Seq:      s.Seq,
		Start:    start,
		End:      end,
		RefID:    s.ID,
		RefSeq:   s,
		RefStart: s.Start,
		RefEnd:   s.End,
	}
	s.SubSeq = append(s.SubSeq, subSeq)
	return subSeq
}

func (s *Seq) CreateFasta(prefix string) {
	fa := osUtil.Create(prefix + ".fa")
	defer simpleUtil.DeferClose(fa)

	fmtUtil.Fprintf(fa, ">%s\n%s\n", s.ID, s.Seq)
}

func main() {
	flag.Parse()
	if *input == "" || *outputDir == "" || *sangerDir == "" {
		flag.PrintDefaults()
		log.Fatal("input is required")
	}

	var inputXlsx = simpleUtil.HandleError(excelize.OpenFile(*input))

	var seqMap = make(map[string]*Seq)
	var seqList []string
	var rows = simpleUtil.HandleError(inputXlsx.GetRows("分段序列"))
	var title []string
	for i, row := range rows {
		if i == 0 {
			title = row
			continue
		}
		var seq = &Seq{}
		for j, cell := range row {
			switch title[j] {
			case "片段名称":
				seq.ID = cell
				seq.RefID = cell[:len(cell)-1]
			case "片段序列":
				seq.Seq = cell
			case "起点":
				seq.Start = stringsUtil.Atoi(cell)
			case "终点":
				seq.End = stringsUtil.Atoi(cell)
			}
		}
		seqMap[seq.ID] = seq
		seqList = append(seqList, seq.ID)
	}

	rows = simpleUtil.HandleError(inputXlsx.GetRows("引物对序列"))
	title = []string{}
	for i, row := range rows {
		if i == 0 {
			title = row
			continue
		}
		var pair = &Seq{}
		var primer1pos [2]int
		var primer2pos [2]int
		for j, cell := range row {
			switch title[j] {
			case "引物对名称":
				pair.ID = cell
				pair.RefID = cell[:len(cell)-1]
			case "引物对序列":
				pair.Seq = cell
			case "有效序列-起点":
				pair.Start = stringsUtil.Atoi(cell)
			case "有效序列-终点":
				pair.End = stringsUtil.Atoi(cell)
			case "左引物-起点":
				primer1pos[0] = stringsUtil.Atoi(cell)
			case "左引物-终点":
				primer1pos[1] = stringsUtil.Atoi(cell)
			case "右引物-起点":
				primer2pos[0] = stringsUtil.Atoi(cell)
			case "右引物-终点":
				primer2pos[1] = stringsUtil.Atoi(cell)
			case "片段-起点":
				pair.RefStart = stringsUtil.Atoi(cell)
			case "片段-终点":
				pair.RefEnd = stringsUtil.Atoi(cell)
			}
		}
		refSeq, ok := seqMap[pair.RefID]
		if !ok {
			log.Fatal("ref not found for:", pair.ID)
		}
		// 检查符合
		if refSeq.Seq[refSeq.Start:refSeq.End][pair.RefStart:pair.RefEnd] != pair.Seq[pair.Start:pair.End] {
			log.Fatal("seq not match for:", pair.ID)
		}
		pair.RefSeq = refSeq
		refSeq.SubSeq = append(refSeq.SubSeq, pair)
		pair.CreateSub(pair.ID+"_1", primer1pos[0], primer1pos[1])
		pair.CreateSub(pair.ID+"_2", primer2pos[0], primer2pos[1])
	}

	simpleUtil.CheckErr(os.MkdirAll(*outputDir, 0755))

	var (
		resultLines  []string
		tracyResults = make(map[string]map[int][2]*tracy.Result)
	)
	for i, id := range seqList {
		seq := seqMap[id]
		slog.Info("seq", "index", i+1, "id", seq.ID, "start", seq.Start, "end", seq.End)

		// 创建 Fasta
		prefix := filepath.Join(*outputDir, id)

		seq.CreateFasta(prefix)
		result := RunTracyBatch(id, prefix, *bin)
		lines := RecordSeq(seq, result, prefix)

		if result != nil {
			tracyResults[id] = result
		}
		resultLines = append(resultLines, lines...)
	}

	// 写入 result
	resultFile := osUtil.Create(filepath.Join(*outputDir, "Result.txt"))
	defer simpleUtil.DeferClose(resultFile)

	fmtUtil.FprintStringArray(resultFile, resultTitle, "\t")
	fmtUtil.FprintStringArray(resultFile, resultLines, "\n")

	// 写入 tracy result
	tracyResultFile := osUtil.Create(filepath.Join(*outputDir, "TracyResult.txt"))
	defer simpleUtil.DeferClose(tracyResultFile)
	fmtUtil.FprintStringArray(tracyResultFile, tracyResultTitle, "\t")
	for id, result := range tracyResults {
		for sangerPairIndex, pairResult := range result {
			for sangerIndex, result := range pairResult {
				fmtUtil.Fprintf(
					tracyResultFile,
					"%s\t%d\t%d\t%s\t%v\t%d\t%d\t%f\n",
					id, sangerPairIndex, sangerIndex,
					result.Status, result.Pass,
					len(result.Variants.Variants),
					result.Variants.HetCount,
					result.AlignResult.BoundMatchRatio,
				)
			}
		}
	}
}
