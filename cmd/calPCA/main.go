package main

import (
	"callAB1/pkg/tracy"
	"flag"
	"fmt"
	"log"
	"log/slog"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"sync"

	"github.com/liserjrqlxue/goUtil/fmtUtil"
	"github.com/liserjrqlxue/goUtil/osUtil"
	"github.com/liserjrqlxue/goUtil/simpleUtil"
	"github.com/liserjrqlxue/goUtil/stringsUtil"
	"github.com/liserjrqlxue/goUtil/textUtil"
	"github.com/xuri/excelize/v2"
)

// flag
var (
	input = flag.String(
		"i",
		"",
		"input 自合.xlsx",
	)
	inputOrder = flag.String(
		"io",
		"",
		"input 引物订购单.xlsx",
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
	renameTxt = flag.String(
		"r",
		"",
		"rename file for CY0130",
	)
	cloneCount = flag.Int(
		"c",
		0,
		"clone count",
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

	if *cloneCount > 0 {
		CloneCountLimit = *cloneCount
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
		if *renameTxt != "" {
			pair.CreateSub(pair.ID, primer1pos[0], primer1pos[1])
		} else {
			pair.CreateSub(pair.ID+"_1", primer1pos[0], primer1pos[1])
			pair.CreateSub(pair.ID+"_2", primer2pos[0], primer2pos[1])
		}
	}

	simpleUtil.CheckErr(os.MkdirAll(*outputDir, 0755))

	var (
		cy0130           bool
		rename           = make(map[string]string)
		resultLines      = make(map[string][][]any)
		tracyStatusLines = make(map[string][][]any)
		seqLines         = make(map[string][]any)
	)
	if *renameTxt != "" {
		cy0130 = true
		rename = simpleUtil.HandleError(textUtil.File2Map(*renameTxt, "\t", false))
	}

	results := make(chan tracyResult, len(seqList)) // 缓冲通道提升性能
	var wg sync.WaitGroup
	for i, id := range seqList {
		wg.Add(1)
		go func(i int, id string) {
			defer wg.Done()
			results <- processSeq(i, id, cy0130, rename, *outputDir, *bin, seqMap)
		}(i, id)
	}

	// 等待所有任务完成并关闭通道
	go func() {
		wg.Wait()
		close(results)
	}()

	for result := range results {
		resultLines[result.id] = result.resultLines
		tracyStatusLines[result.id] = result.statusLines
		seqLines[result.id] = result.seqLines
	}

	// 写入 result
	resultFormat := "%s\t%s\t%3d-%3d\t%3d-%3d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f%%\t%.4f%%\t%.4f%%\t%.4f%%\t%4.f%%\t%.4f%%\t%.4f%%\t%.4f%%\n"
	WriteSlice(filepath.Join(*outputDir, "Result.txt"), resultFormat, ResultTitle, seqList, resultLines)

	// 写入 tracy result
	tracyFormat := "%s\t%s\t%d\t%s\t%v\t%d\t%d\t%f\n"
	WriteSlice(filepath.Join(*outputDir, "TracyResult.txt"), tracyFormat, tracyStatusTitle, seqList, tracyStatusLines)

	// 写入 excel
	simpleUtil.HandleError(inputXlsx.NewSheet("Sanger结果"))
	inputXlsx.SetSheetRow("Sanger结果", "A1", &ResultTitle)
	primerACC := make(map[string][2]float64)
	row := 2
	for _, id := range seqList {
		for _, line := range resultLines[id] {
			pID := line[1].(string)
			acc := line[19].(float64)
			geomMeanAcc := line[21].(float64)
			primerACC[pID] = [2]float64{acc, geomMeanAcc}
			inputXlsx.SetSheetRow("Sanger结果", fmt.Sprintf("A%d", row), &line)
			row++
		}
	}

	simpleUtil.HandleError(inputXlsx.NewSheet("片段结果"))
	inputXlsx.SetSheetRow("片段结果", "A1", &SeqTitle)
	row = 2
	for _, id := range seqList {
		line := seqLines[id]
		inputXlsx.SetSheetRow("片段结果", fmt.Sprintf("A%d", row), &line)
		row++
	}

	if *inputOrder != "" {
		inputOderXlsx := simpleUtil.HandleError(excelize.OpenFile(*inputOrder))
		rows := simpleUtil.HandleError(inputOderXlsx.GetRows("拼接引物板"))
		simpleUtil.HandleError(inputXlsx.NewSheet("拼接引物板"))
		for i, row := range rows {
			for j, cell := range row {
				inputXlsx.SetCellValue(
					"拼接引物板",
					simpleUtil.HandleError(excelize.CoordinatesToCellName(j+1, i+1)),
					cell,
				)
				pID := strings.Split(cell, "\n")[0]
				accs, ok := primerACC[pID]
				if ok {
					inputXlsx.SetCellValue(
						"拼接引物板",
						simpleUtil.HandleError(excelize.CoordinatesToCellName(j+1, i+1+len(rows)+1)),
						accs[0],
					)
					inputXlsx.SetCellValue(
						"拼接引物板",
						simpleUtil.HandleError(excelize.CoordinatesToCellName(j+1, i+1+len(rows)*2+2)),
						accs[1],
					)
				} else {
					inputXlsx.SetCellValue(
						"拼接引物板",
						simpleUtil.HandleError(excelize.CoordinatesToCellName(j+1, i+1+len(rows)+1)),
						cell,
					)
					inputXlsx.SetCellValue(
						"拼接引物板",
						simpleUtil.HandleError(excelize.CoordinatesToCellName(j+1, i+1+len(rows)*2+2)),
						cell,
					)
				}
			}
			inputXlsx.SetCellValue(
				"拼接引物板",
				simpleUtil.HandleError(excelize.CoordinatesToCellName(1, 1)),
				"板位",
			)
			inputXlsx.SetCellValue(
				"拼接引物板",
				simpleUtil.HandleError(excelize.CoordinatesToCellName(1, 1+len(rows)+1)),
				"平均准确率(%)",
			)
			inputXlsx.SetCellValue(
				"拼接引物板",
				simpleUtil.HandleError(excelize.CoordinatesToCellName(1, 1+len(rows)*2+2)),
				"参考单步准确率(%)",
			)
		}
	}

	inputXlsx.SaveAs(*outputDir + ".Sanger结果.xlsx")

}

type tracyResult struct {
	id          string
	statusLines [][]any
	resultLines [][]any
	seqLines    []any
}
type Variant struct {
	GeneID    string
	CloneID   string
	SangerID  string
	VariantID string
	Pos       int
	Variant   *tracy.Variant
}

var VariantStringTitle = append(
	[]string{
		"CloneID",
		"SangerID",
		"VariantID",
		"Pos",
	},
	tracy.VariantTitle...,
)

func (v *Variant) String() string {
	return fmt.Sprintf(
		"%s\t%s\t%s\t%d\t%s",
		v.CloneID, v.SangerID, v.VariantID,
		v.Pos,
		v.Variant,
	)
}

type CloneVariant struct {
	GeneID      string
	CloneID     string
	VariantID   string
	Pos         int
	SangerID    []string
	SangerCount int
	Variant     *tracy.Variant
}

var CloneVariantStringTitle = []string{
	"GeneID",
	"CloneID",
	"VariantID",
	"Pos",
	"SangerCount",
	"Chr",
	"Pos",
	"Ref",
	"Alt",
	"Type",
}

func (v *CloneVariant) String() string {
	return fmt.Sprintf(
		"%s\t%s\t%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s",
		v.GeneID, v.CloneID, v.VariantID,
		v.Pos, v.SangerCount,
		v.Variant.Chr,
		v.Variant.Pos,
		v.Variant.Ref,
		v.Variant.Alt,
		v.Variant.Type,
	)
}

type VariantSet struct {
	GeneID     string
	VariantID  string
	Pos        int
	CloneID    []string
	CloneCount int
	Variant    *tracy.Variant
}

var VariantSetStringTitle = []string{
	"GeneID",
	"VariantID",
	"Pos",
	"CloneCount",
	"Chr",
	"Pos",
	"Ref",
	"Alt",
	"Type",
}

func (v *VariantSet) String() string {
	return fmt.Sprintf(
		"%s\t%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s",
		v.GeneID, v.VariantID,
		v.Pos, v.CloneCount,
		v.Variant.Chr,
		v.Variant.Pos,
		v.Variant.Ref,
		v.Variant.Alt,
		v.Variant.Type,
	)
}

type PosVariantSet struct {
	GeneID       string
	Pos          int
	VariantID    map[string]int
	variantCount int
	Variant      *tracy.Variant
}

var PosVariantSetStringTitle = []string{
	"GeneID",
	"Pos",
	"VariantCount",
	"Chr",
	"Pos",
}

func (v *PosVariantSet) String() string {
	return fmt.Sprintf(
		"%s\t%d\t%d\t%s\t%d",
		v.GeneID,
		v.Pos, v.variantCount,
		v.Variant.Chr,
		v.Variant.Pos,
	)
}

func processSeq(i int, id string, cy0130 bool, rename map[string]string, outputDir string, bin string, seqMap map[string]*Seq) tracyResult {
	seq := seqMap[id]
	slog.Info("seq", "index", i+1, "id", seq.ID, "start", seq.Start, "end", seq.End)

	prefix := filepath.Join(outputDir, id)
	seq.CreateFasta(prefix)

	var result map[string][2]*tracy.Result
	if cy0130 {
		result = RunTracyBatchCy0130(rename[id], prefix, bin)
	} else {
		result = RunTracyBatch(id, prefix, bin)
	}

	var variants []*Variant
	for cloneID, results := range result {
		for i, result := range results {
			if result == nil || result.Variants == nil {
				continue
			}
			for _, variant := range result.Variants.Variants {
				variants = append(variants, &Variant{
					GeneID:    id,
					CloneID:   cloneID,
					SangerID:  fmt.Sprintf("%s_%d", cloneID, i+1),
					VariantID: fmt.Sprintf("%s_%d_%s_%s", variant.Chr, variant.Pos, variant.Ref, variant.Alt),
					Pos:       variant.Pos,
					Variant:   variant,
				})
			}
		}
	}
	// 按Pos排序
	sort.Slice(
		variants,
		func(i, j int) bool {
			if variants[i].Pos == variants[j].Pos {
				if variants[i].VariantID == variants[j].VariantID {
					return variants[i].SangerID < variants[j].SangerID
				}
				return variants[i].VariantID < variants[j].VariantID
			}
			return variants[i].Pos < variants[j].Pos
		},
	)
	out := osUtil.Create(filepath.Join(outputDir, id+".variant.raw.txt"))
	fmtUtil.FprintStringArray(out, VariantStringTitle, "\t")
	// 合并同一克隆的变异
	var cloneVariants = make(map[string]*CloneVariant)
	for _, variant := range variants {
		fmtUtil.Fprintln(out, variant)
		key := variant.CloneID + "_" + variant.VariantID
		cv, ok := cloneVariants[key]
		if ok {
			cv.SangerID = append(cv.SangerID, variant.SangerID)
			cv.SangerCount = len(cv.SangerID)
		} else {
			cv = &CloneVariant{
				GeneID:      variant.GeneID,
				CloneID:     variant.CloneID,
				VariantID:   variant.VariantID,
				Pos:         variant.Pos,
				SangerID:    []string{variant.SangerID},
				SangerCount: 1,
				Variant:     variant.Variant,
			}
		}
		cloneVariants[key] = cv
	}
	simpleUtil.CheckErr(out.Close())

	// 合并同一变异
	var variantsSet = make(map[string]*VariantSet)
	var cloneVariantsArray []*CloneVariant
	for _, cl := range cloneVariants {
		cloneVariantsArray = append(cloneVariantsArray, cl)
		key := cl.VariantID
		vs, ok := variantsSet[key]
		if ok {
			vs.CloneID = append(vs.CloneID, cl.CloneID)
			vs.CloneCount = len(vs.CloneID)
		} else {
			vs = &VariantSet{
				GeneID:     cl.GeneID,
				VariantID:  cl.VariantID,
				Pos:        cl.Pos,
				CloneID:    []string{cl.CloneID},
				CloneCount: 1,
				Variant:    cl.Variant,
			}
		}
		variantsSet[key] = vs
	}

	// 按Pos排序
	sort.Slice(
		cloneVariantsArray,
		func(i, j int) bool {
			if cloneVariantsArray[i].Pos == cloneVariantsArray[j].Pos {
				if cloneVariantsArray[i].VariantID == cloneVariantsArray[j].VariantID {
					return cloneVariantsArray[i].CloneID < cloneVariantsArray[j].CloneID
				}
				return cloneVariantsArray[i].VariantID < cloneVariantsArray[j].VariantID
			}
			return cloneVariantsArray[i].Pos < cloneVariantsArray[j].Pos
		},
	)
	out = osUtil.Create(filepath.Join(outputDir, id+".variant.clone.txt"))
	fmtUtil.FprintStringArray(out, CloneVariantStringTitle, "\t")
	for _, cv := range cloneVariantsArray {
		fmtUtil.Fprintln(out, cv)
	}
	simpleUtil.CheckErr(out.Close())

	// 合并同一位置
	var posVariantsSet = make(map[int]*PosVariantSet)
	var variantsSetArray []*VariantSet
	for _, vs := range variantsSet {
		variantsSetArray = append(variantsSetArray, vs)
		key := vs.Pos
		pvs, ok := posVariantsSet[key]
		if ok {
			pvs.VariantID[vs.VariantID]++
			pvs.variantCount++
		} else {
			pvs = &PosVariantSet{
				GeneID:       vs.GeneID,
				Pos:          vs.Pos,
				VariantID:    map[string]int{vs.VariantID: 1},
				variantCount: 1,
				Variant:      vs.Variant,
			}
		}
		posVariantsSet[key] = pvs
	}
	// 按Pos排序
	sort.Slice(
		variantsSetArray,
		func(i, j int) bool {
			if variantsSetArray[i].Pos == variantsSetArray[j].Pos {
				return variantsSetArray[i].VariantID < variantsSetArray[j].VariantID
			}
			return variantsSetArray[i].Pos < variantsSetArray[j].Pos
		},
	)
	out = osUtil.Create(filepath.Join(outputDir, id+".variant.set.txt"))
	fmtUtil.FprintStringArray(out, VariantSetStringTitle, "\t")
	for _, vs := range variantsSetArray {
		fmtUtil.Fprintln(out, vs)
	}
	simpleUtil.CheckErr(out.Close())

	var posSet []int
	for pos := range posVariantsSet {
		posSet = append(posSet, pos)
	}
	sort.Ints(posSet)
	out = osUtil.Create(filepath.Join(outputDir, id+".variant.pos.txt"))
	fmtUtil.FprintStringArray(out, PosVariantSetStringTitle, "\t")
	for _, pos := range posSet {
		psv := posVariantsSet[pos]
		fmtUtil.Fprintln(out, psv)
	}
	simpleUtil.CheckErr(out.Close())

	return tracyResult{
		id:          id,
		statusLines: GetTracyStatusLines(id, result),
		resultLines: RecordSeqPrimer(seq, result, prefix),
		seqLines:    RecordSeq(seq, result, prefix),
	}
}
