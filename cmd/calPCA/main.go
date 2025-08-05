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
	"github.com/liserjrqlxue/goUtil/textUtil"
	"github.com/liserjrqlxue/version"
	"github.com/samber/lo"
	"github.com/xuri/excelize/v2"
)

// os
var (
	ex, _  = os.Executable()
	exPath = filepath.Dir(ex)
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
	maxQual = flag.Int(
		"q",
		0,
		"filter qual < MaxQual",
	)
	override = flag.Bool(
		"w",
		false,
		"override tracy result",
	)
	fix = flag.Bool(
		"fix",
		false,
		"run chemical supplements",
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

	rIDs []string

	CloneHit int // 正确克隆个数
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

	MaxQual = *maxQual

	var (
		cy0130            bool
		rename            = make(map[string]string)
		resultLines       = make(map[string][][]any)
		tracyStatusLines  = make(map[string][][]any)
		segmentLines      = make(map[string][]any) // 片段结果
		cloneVariantLines = make(map[string][][]any)
		setVariantLines   = make(map[string][][]any)

		xlsx  = simpleUtil.HandleError(excelize.OpenFile(*input))
		sheet string

		geneMap, geneList       = LoadRawSequence(xlsx, "原始序列")
		segmentMap, segmentList = LoadSegmentSequence(xlsx, "分段序列", geneMap)

		sangerResult = *outputDir + ".Sanger结果.xlsx"
		// batchName    = strings.TrimSuffix(filepath.Base(*input), "自合.xlsx")
		// csPrefix     = filepath.Join(filepath.Dir(*outputDir), batchName)
		template = filepath.Join(exPath, "AZENTA.xlsx")
	)

	LoadPrimerPairSequence(xlsx, "引物对序列", segmentMap)

	simpleUtil.CheckErr(os.MkdirAll(*outputDir, 0755))

	if *renameTxt != "" {
		cy0130 = true
		if osUtil.FileExists(*renameTxt) {
			rename = simpleUtil.HandleError(textUtil.File2Map(*renameTxt, "\t", false))
		} else {
			for _, id := range segmentList {
				// rename[id] = strings.Replace(id, "A", "-A", 1)
				rename[id] = id
			}
		}
		slog.Info("RENAME", "rename", rename)
	}

	results := make(chan tracyResult, len(segmentList)) // 缓冲通道提升性能
	var wg sync.WaitGroup
	for i, id := range segmentList {
		renameID, ok := rename[id]
		if !ok {
			slog.Error("Skip", "i", i, "id", id)
		} else {
			wg.Add(1)
			go func(i int, id string) {
				defer wg.Done()
				results <- processSeq(i, id, cy0130, renameID, *outputDir, *bin, segmentMap, *override)
			}(i, id)
		}
	}

	// 等待所有任务完成并关闭通道
	go func() {
		wg.Wait()
		close(results)
	}()

	var (
		ratios     [3][]float64
		batchVerif = &Verification{Status: "合格"}
	)
	for result := range results {
		resultLines[result.ID] = result.resultLines
		tracyStatusLines[result.ID] = result.statusLines
		segmentLines[result.ID] = result.segmentLines
		cloneVariantLines[result.ID] = result.cloneVariantLines
		setVariantLines[result.ID] = result.setVariantLines

		ratios[0] = append(ratios[0], result.SnvRatio)
		ratios[1] = append(ratios[1], result.InsRatio)
		ratios[2] = append(ratios[2], result.DelRatio)
	}

	// 写入 result
	resultFormat := "%s\t%s\t%3d-%3d\t%3d-%3d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f%%\t%.4f%%\t%.4f%%\t%.4f%%\t%4.f%%\t%.4f%%\t%.4f%%\t%.4f%%\n"
	WriteSlice(filepath.Join(*outputDir, "Result.txt"), resultFormat, ResultTitle, segmentList, resultLines)

	// 写入 tracy result
	tracyFormat := "%s\t%s\t%d\t%s\t%v\t%d\t%d\t%f\n"
	WriteSlice(filepath.Join(*outputDir, "TracyResult.txt"), tracyFormat, tracyStatusTitle, segmentList, tracyStatusLines)
	WriteSliceSheet(xlsx, "Sanger统计", tracyStatusTitle, segmentList, tracyStatusLines)

	// 写入 excel
	simpleUtil.HandleError(xlsx.NewSheet("Sanger结果"))
	xlsx.SetSheetRow("Sanger结果", "A1", &ResultTitle)
	primerACC := make(map[string][2]float64)
	row := 2
	for _, id := range segmentList {
		for _, line := range resultLines[id] {
			pID := line[1].(string)
			valid := line[7].(int)
			acc := line[19].(float64)
			geomMeanAcc := line[21].(float64)

			if valid > 0 {
				primerACC[pID] = [2]float64{acc, geomMeanAcc}
			}
			xlsx.SetSheetRow("Sanger结果", fmt.Sprintf("A%d", row), &line)
			row++
		}
	}

	simpleUtil.HandleError(xlsx.NewSheet("片段结果"))
	xlsx.SetSheetRow("片段结果", "A1", &SeqTitle)
	row = 2
	for _, id := range segmentList {
		line := segmentLines[id]
		xlsx.SetSheetRow("片段结果", fmt.Sprintf("A%d", row), &line)
		row++
	}

	AddSequencingResultPlate(xlsx, "测序结果板位图", geneList, geneMap)

	simpleUtil.HandleError(xlsx.NewSheet("Clone变异结果"))
	xlsx.SetSheetRow("Clone变异结果", "A1", &CloneVariantTitle)
	row = 2
	for _, id := range segmentList {
		lines := cloneVariantLines[id]
		for _, line := range lines {
			xlsx.SetSheetRow("Clone变异结果", fmt.Sprintf("A%d", row), &line)
			row++
		}
	}

	simpleUtil.HandleError(xlsx.NewSheet("变异统计"))
	xlsx.SetSheetRow("变异统计", "A1", &SetVariantTitle)
	row = 2
	for _, id := range segmentList {
		lines := setVariantLines[id]
		for _, line := range lines {
			xlsx.SetSheetRow("变异统计", fmt.Sprintf("A%d", row), &line)
			row++
		}
	}

	// 拼接引物板
	if *inputOrder != "" {
		sheet = "拼接引物板"
		addSplicedPrimerPlate(xlsx, sheet, *inputOrder, primerACC)
	}

	sheet = "批次统计"
	batchVerif.BatchValidation(ratios[0], ratios[1], ratios[2])
	addBatchStats(xlsx, sheet, batchVerif)

	log.Printf("SaveAs(%s)", sangerResult)
	simpleUtil.CheckErr(xlsx.SaveAs(sangerResult))

	// 化学补充
	if *fix {
		runFix(*outputDir, template, segmentList, segmentMap)
	}
}

type tracyResult struct {
	ID  string
	Seq *Seq

	statusLines       [][]any
	resultLines       [][]any
	segmentLines      []any
	cloneVariantLines [][]any
	setVariantLines   [][]any

	// (%)
	SnvRatio float64
	InsRatio float64
	DelRatio float64
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

func processSeq(i int, id string, cy0130 bool, renameID, outputDir string, bin string, seqMap map[string]*Seq, override bool) tracyResult {
	seq := seqMap[id]
	slog.Info("seq", "index", i+1, "id", seq.ID, "start", seq.Start, "end", seq.End)

	prefix := filepath.Join(outputDir, id)
	seq.CreateFasta(prefix)

	var result map[string][2]*tracy.Result
	if cy0130 {
		result = RunTracyBatchCy0130(renameID, prefix, bin, override)
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

	seqLines, vTypePercent := RecordSeq(seq, result, prefix)

	return tracyResult{
		ID:  id,
		Seq: seq,

		statusLines:       GetTracyStatusLines(id, result),
		resultLines:       RecordSeqPrimer(seq, result, prefix),
		segmentLines:      seqLines,
		cloneVariantLines: cloneVariantsLines,
		setVariantLines:   setVariantLines,

		SnvRatio: vTypePercent["SNV"],
		InsRatio: vTypePercent["Insertion"],
		DelRatio: vTypePercent["Deletion"],
	}
}
