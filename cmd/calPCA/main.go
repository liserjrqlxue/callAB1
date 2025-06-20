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
	"strconv"
	"strings"
	"sync"

	"github.com/liserjrqlxue/goUtil/fmtUtil"
	"github.com/liserjrqlxue/goUtil/osUtil"
	"github.com/liserjrqlxue/goUtil/simpleUtil"
	"github.com/liserjrqlxue/goUtil/stringsUtil"
	"github.com/liserjrqlxue/goUtil/textUtil"
	"github.com/liserjrqlxue/version"
	"github.com/samber/lo"
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
	version.LogVersion()
	flag.Parse()
	if *outputDir == "" {
		flag.PrintDefaults()
		log.Fatal("-o is required")
	}
	if *input == "" {
		*input = filepath.Join(*outputDir, "自合.xlsx")
	}
	if *inputOrder == "" {
		*inputOrder = filepath.Join(*outputDir, "引物订购单_BOM.xlsx")
	}
	if *renameTxt == "" {
		*renameTxt = filepath.Join(*outputDir, "rename.txt")
	}
	if *sangerDir == "" {
		*renameTxt = filepath.Join(*outputDir, "ab1")
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
				if cell == "" {
					log.Fatalf("分段序列 片段名称(%d,%d) 为空", i+1, j+1)
				}
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
				pair.RefID = strings.TrimSuffix(cell[:len(cell)-1], "_")
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
			log.Fatalf("ref not found for:[id:%s,ref:%s],[keys:%+v]", pair.ID, pair.RefID, lo.Keys(seqMap))
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
		cy0130            bool
		rename            = make(map[string]string)
		resultLines       = make(map[string][][]any)
		tracyStatusLines  = make(map[string][][]any)
		seqLines          = make(map[string][]any)
		cloneVariantLines = make(map[string][][]any)
		setVariantLines   = make(map[string][][]any)
	)
	if *renameTxt != "" {
		cy0130 = true
		if osUtil.FileExists(*renameTxt) {
			rename = simpleUtil.HandleError(textUtil.File2Map(*renameTxt, "\t", false))
		} else {
			for _, id := range seqList {
				// rename[id] = strings.Replace(id, "A", "-A", 1)
				rename[id] = id
			}
		}
		slog.Info("RENAME", "rename", rename)
	}

	results := make(chan tracyResult, len(seqList)) // 缓冲通道提升性能
	var wg sync.WaitGroup
	for i, id := range seqList {
		renameID, ok := rename[id]
		if !ok {
			slog.Error("Skip", "i", i, "id", id)
		} else {
			wg.Add(1)
			go func(i int, id string) {
				defer wg.Done()
				results <- processSeq(i, id, cy0130, renameID, *outputDir, *bin, seqMap)
			}(i, id)
		}
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
		cloneVariantLines[result.id] = result.cloneVariantLines
		setVariantLines[result.id] = result.setVariantLines
	}

	// 写入 result
	resultFormat := "%s\t%s\t%3d-%3d\t%3d-%3d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f%%\t%.4f%%\t%.4f%%\t%.4f%%\t%4.f%%\t%.4f%%\t%.4f%%\t%.4f%%\n"
	WriteSlice(filepath.Join(*outputDir, "Result.txt"), resultFormat, ResultTitle, seqList, resultLines)

	// 写入 tracy result
	tracyFormat := "%s\t%s\t%d\t%s\t%v\t%d\t%d\t%f\n"
	WriteSlice(filepath.Join(*outputDir, "TracyResult.txt"), tracyFormat, tracyStatusTitle, seqList, tracyStatusLines)
	WriteSliceSheet(inputXlsx, "Sanger统计", tracyStatusTitle, seqList, tracyStatusLines)

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

	simpleUtil.HandleError(inputXlsx.NewSheet("Clone变异结果"))
	inputXlsx.SetSheetRow("Clone变异结果", "A1", &CloneVariantTitle)
	row = 2
	for _, id := range seqList {
		lines := cloneVariantLines[id]
		for _, line := range lines {
			inputXlsx.SetSheetRow("Clone变异结果", fmt.Sprintf("A%d", row), &line)
			row++
		}
	}

	simpleUtil.HandleError(inputXlsx.NewSheet("变异统计"))
	inputXlsx.SetSheetRow("变异统计", "A1", &SetVariantTitle)
	row = 2
	for _, id := range seqList {
		lines := setVariantLines[id]
		for _, line := range lines {
			inputXlsx.SetSheetRow("变异统计", fmt.Sprintf("A%d", row), &line)
			row++
		}
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
	id                string
	statusLines       [][]any
	resultLines       [][]any
	seqLines          []any
	cloneVariantLines [][]any
	setVariantLines   [][]any
}
type Variant struct {
	GeneID       string
	CloneID      string
	SangerID     string
	SangerStatus string
	SangerPass   bool
	VariantID    string
	Pos          int
	Variant      *tracy.Variant
}

var VariantStringTitle = append(
	tracy.VariantTitle,
	"CloneID",
	"SangerID",
	"SangerStatus",
	"SangerPass",
)

func (v *Variant) String() string {
	return fmt.Sprintf(
		"%s\t%s\t%s\t%s\t%t",
		v.Variant,
		v.CloneID,
		v.SangerID,
		v.SangerStatus,
		v.SangerPass,
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

	Qual     []int
	Filter   []string
	Genotype []string
}

func (v *CloneVariant) String() string {
	return fmt.Sprintf(
		"%s\t%d\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%s",
		v.Variant.Chr,
		v.Variant.Pos,
		v.Variant.Ref,
		v.Variant.Alt,
		v.Variant.Type,
		v.CloneID,
		v.SangerCount,
		lo.Max(v.Qual),
		strings.Join(lo.Uniq(v.Filter), ","),
		strings.Join(lo.Uniq(v.Genotype), ","),
	)
}

type VariantSet struct {
	GeneID     string
	VariantID  string
	Pos        int
	CloneID    []string
	CloneCount int
	ClonePass  int
	Variant    *tracy.Variant

	Qual     []int
	Filter   []string
	Genotype []string
}

func (v *VariantSet) String() string {
	return fmt.Sprintf(
		"%s\t%d\t%s\t%s\t%s\t%d\t%d\t%f\t%d\t%s\t%s",
		v.Variant.Chr,
		v.Variant.Pos,
		v.Variant.Ref,
		v.Variant.Alt,
		v.Variant.Type,
		v.CloneCount,
		v.ClonePass,
		float64(v.CloneCount)/float64(v.ClonePass),
		lo.Max(v.Qual),
		strings.Join(lo.Uniq(v.Filter), ","),
		strings.Join(lo.Uniq(v.Genotype), ","),
	)
}

type PosVariantSet struct {
	GeneID       string
	Pos          int
	VariantID    map[string]int
	VariantCount int
	ClonePass    int
	Variant      *tracy.Variant
}

var PosVariantSetStringTitle = []string{
	"Chr",
	"Pos",
	"VariantCount",
	"ClonePass",
	"PosRatio",
}

func (v *PosVariantSet) String() string {
	return fmt.Sprintf(
		"%s\t%d\t%d\t%d\t%f",
		v.Variant.Chr,
		v.Variant.Pos,
		v.VariantCount,
		v.ClonePass,
		float64(v.VariantCount)/float64(v.ClonePass),
	)
}

func processSeq(i int, id string, cy0130 bool, renameID, outputDir string, bin string, seqMap map[string]*Seq) tracyResult {
	seq := seqMap[id]
	slog.Info("seq", "index", i+1, "id", seq.ID, "start", seq.Start, "end", seq.End)

	prefix := filepath.Join(outputDir, id)
	seq.CreateFasta(prefix)

	var result map[string][2]*tracy.Result
	if cy0130 {
		result = RunTracyBatchCy0130(renameID, prefix, bin)
	} else {
		result = RunTracyBatch(id, prefix, bin)
	}

	var (
		variants           []*Variant
		cloneVariantsLines [][]any
		setVariantLines    [][]any
	)
	var clonePass int // 有效克隆数
	for cloneID, results := range result {
		pass := false
		for i, result := range results {
			if result == nil || result.Variants == nil {
				continue
			}
			if result.Pass {
				pass = true
			}
			for _, variant := range result.Variants.Variants {
				variants = append(variants, &Variant{
					GeneID:       id,
					CloneID:      cloneID,
					SangerID:     strconv.Itoa(i + 1),
					SangerStatus: result.Status,
					SangerPass:   result.Pass,
					VariantID:    fmt.Sprintf("%s_%d_%s_%s", variant.Chr, variant.Pos, variant.Ref, variant.Alt),
					Pos:          variant.Pos,
					Variant:      variant,
				})
			}
		}
		if pass {
			clonePass++
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

		if !variant.SangerPass {
			continue
		}

		key := variant.CloneID + "_" + variant.VariantID
		cv, ok := cloneVariants[key]
		if ok {
			cv.SangerID = append(cv.SangerID, variant.SangerID)
			cv.Qual = append(cv.Qual, variant.Variant.Qual)
			cv.Filter = append(cv.Filter, variant.Variant.Filter)
			cv.Genotype = append(cv.Genotype, variant.Variant.Type)
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
				Qual:        []int{variant.Variant.Qual},
				Filter:      []string{variant.Variant.Filter},
				Genotype:    []string{variant.Variant.Genotype},
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
			vs.Qual = append(vs.Qual, cl.Qual...)
			vs.Filter = append(vs.Filter, cl.Filter...)
			vs.Genotype = append(vs.Genotype, vs.Genotype...)
			vs.CloneCount = len(vs.CloneID)
		} else {
			vs = &VariantSet{
				GeneID:     cl.GeneID,
				VariantID:  cl.VariantID,
				Pos:        cl.Pos,
				CloneID:    []string{cl.CloneID},
				CloneCount: 1,
				ClonePass:  clonePass,
				Variant:    cl.Variant,
				Qual:       cl.Qual,
				Filter:     cl.Filter,
				Genotype:   cl.Genotype,
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
	fmtUtil.FprintStringArray(out, CloneVariantTitle, "\t")
	for _, cv := range cloneVariantsArray {
		fmtUtil.Fprintln(out, cv)
		cvLine := []any{
			cv.Variant.Chr,
			cv.Variant.Pos,
			cv.Variant.Ref,
			cv.Variant.Alt,
			cv.Variant.Type,
			cv.CloneID,
			cv.SangerCount,
			lo.Max(cv.Qual),
			strings.Join(lo.Uniq(cv.Filter), ","),
			strings.Join(lo.Uniq(cv.Genotype), ","),
		}
		cloneVariantsLines = append(cloneVariantsLines, cvLine)
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
			pvs.VariantID[vs.VariantID] += vs.CloneCount
			pvs.VariantCount += vs.CloneCount
		} else {
			pvs = &PosVariantSet{
				GeneID:       vs.GeneID,
				Pos:          vs.Pos,
				VariantID:    map[string]int{vs.VariantID: vs.CloneCount},
				VariantCount: vs.CloneCount,
				ClonePass:    clonePass,
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
	fmtUtil.FprintStringArray(out, SetVariantTitle, "\t")
	for _, vs := range variantsSetArray {
		fmtUtil.Fprintln(out, vs)
		vsLine := []any{
			vs.Variant.Chr,
			vs.Variant.Pos,
			vs.Variant.Ref,
			vs.Variant.Alt,
			vs.Variant.Type,
			vs.CloneCount,
			vs.ClonePass,
			float64(vs.CloneCount) / float64(vs.ClonePass),
			lo.Max(vs.Qual),
			strings.Join(lo.Uniq(vs.Filter), ","),
			strings.Join(lo.Uniq(vs.Genotype), ","),
		}
		setVariantLines = append(setVariantLines, vsLine)
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
		id:                id,
		statusLines:       GetTracyStatusLines(id, result),
		resultLines:       RecordSeqPrimer(seq, result, prefix),
		seqLines:          RecordSeq(seq, result, prefix),
		cloneVariantLines: cloneVariantsLines,
		setVariantLines:   setVariantLines,
	}
}
