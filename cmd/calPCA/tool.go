package main

import (
	"callAB1/pkg/tracy"
	"fmt"
	"log"
	"log/slog"
	"math"
	"os"
	"path/filepath"
	"regexp"
	"sort"
	"strconv"
	"strings"

	"github.com/liserjrqlxue/goUtil/fmtUtil"
	"github.com/liserjrqlxue/goUtil/osUtil"
	"github.com/liserjrqlxue/goUtil/simpleUtil"
	"github.com/liserjrqlxue/goUtil/stringsUtil"
	"github.com/samber/lo"
	"github.com/xuri/excelize/v2"
)

func RunTracy(id, tag, prefix, bin, sangerIndex string, result map[string][2]*tracy.Result) {
	sangerPrefix := fmt.Sprintf("%s_%s", prefix, sangerIndex)
	path1 := filepath.Join(*sangerDir, fmt.Sprintf("BGA-%s-TY1-%s-T7.ab1", tag, sangerIndex))
	path2 := filepath.Join(*sangerDir, fmt.Sprintf("BGA-%s-TY1-%s-T7-Term.ab1", tag, sangerIndex))

	ok1 := osUtil.FileExists(path1)
	ok2 := osUtil.FileExists(path2)
	if ok1 && ok2 {
		// run tracy
		result1, result2, err := tracy.RunPair(bin, prefix+".fa", path1, path2, sangerPrefix)
		simpleUtil.CheckErr(err)
		result[sangerIndex] = [2]*tracy.Result{&result1, &result2}
	} else {
		slog.Debug("sanger pair no found", "id", id, "sangerIndex", sangerIndex, "path1", path1, "ok1", ok1, "path2", path2, "ok2", ok2)
		if ok1 || ok2 {
			log.Fatalf("sanger pair fail:[%s:%v,%s:%v]", path1, ok1, path2, ok2)
		}
	}
}

func RunTracyCY0130(id, prefix, bin, path string, result map[string][2]*tracy.Result) {
	sangerIndex := strings.TrimSuffix(filepath.Base(path), ".ab1")
	sangerPrefix := fmt.Sprintf("%s_%s", prefix, sangerIndex)

	ok := osUtil.FileExists(path)
	if !ok {
		slog.Info("sanger file no found", "id", id, "sangerIndex", sangerIndex, "path", path, "ok", ok)
		return
	}

	// run tracy
	result1, err := tracy.RunSingle(bin, prefix+".fa", path, sangerPrefix, false)
	if err != nil {
		slog.Error("RunSingle", "id", id, "sangerIndex", sangerIndex, "err", err)
	}
	result[sangerIndex] = [2]*tracy.Result{&result1}
}

// 遍历分析sanger文件 -> result
func RunTracyBatch(id, prefix, bin string) map[string][2]*tracy.Result {
	result := make(map[string][2]*tracy.Result)
	tag := []byte(strings.Split(id, "_")[1][:4])
	tag[2] = 'B'
	for sangerIndex := 1; sangerIndex <= CloneCountLimit; sangerIndex++ {
		RunTracy(id, string(tag), prefix, bin, strconv.Itoa(sangerIndex), result)
	}
	return result
}

// 遍历分析sanger文件 -> result
func RunTracyBatchCy0130(id, prefix, bin string) map[string][2]*tracy.Result {
	result := make(map[string][2]*tracy.Result)

	files := simpleUtil.HandleError(
		filepath.Glob(fmt.Sprintf("%s/%s*.ab1", *sangerDir, id)),
	)
	files2 := simpleUtil.HandleError(
		filepath.Glob(fmt.Sprintf("%s/*/%s*.ab1", *sangerDir, id)),
	)
	files = append(files, files2...)

	for _, path := range files {
		RunTracyCY0130(id, prefix, bin, path, result)
	}

	return result
}

func RecordVariant(v *tracy.Variant, out *os.File, index, id, sangerPairIndex string, sangerIndex, start, end, hetCount int, boundMatchRatio float64, variantSet map[string]bool, pass bool) {
	changeID := strings.ReplaceAll(sangerPairIndex, "-", "_")
	if v.Pos > start && v.Pos <= end {
		fmtUtil.Fprintf(
			out,
			"%s\t%s\t%3d-%3d\t%s.%d\t%f\t%d\t%v\t%s\n",
			index,
			id,
			start, end,
			sangerPairIndex, sangerIndex,
			boundMatchRatio, hetCount,
			pass,
			v,
		)
		if pass {
			if v.Qual >= MaxQual {
				key := fmt.Sprintf("%s-%d-%s-%s-%s", changeID, v.Pos, v.Ref, v.Alt, v.Type)
				variantSet[key] = true
			}
		}
	}
}

func Record1Result(primer *Seq, result map[string][2]*tracy.Result, out *os.File, index, sangerIndex string, offset int, variantSet map[string]bool) (keep bool) {
	var (
		id    = primer.ID
		start = primer.Start + offset
		end   = primer.End + offset

		result1 = result[sangerIndex][0]
		result2 = result[sangerIndex][1]
	)

	if (result1 != nil && result1.Pass) || (result2 != nil && result2.Pass) {
		keep = true
	}

	if result1 != nil && result1.Variants != nil {
		for _, v := range result1.Variants.Variants {
			if v.SV {
				v.Type = "SV"
			}
			RecordVariant(
				v, out, index, id, sangerIndex, 1, start, end,
				result1.Variants.HetCount,
				result1.AlignResult.BoundMatchRatio,
				variantSet,
				result1.Pass,
			)
		}
	}
	if result2 != nil && result2.Variants != nil {
		for _, v := range result2.Variants.Variants {
			RecordVariant(
				v, out, index, id, sangerIndex, 2, start, end,
				result2.Variants.HetCount,
				result2.AlignResult.BoundMatchRatio,
				variantSet,
				result2.Pass,
			)
		}
	}
	return
}

func RecordPrimer(primer *Seq, result map[string][2]*tracy.Result, out *os.File, index string, offset int) (resultLine []any) {
	var (
		length            = primer.End - primer.Start // 长度
		n                 = 0                         // sanger 个数
		ln                = 0                         // length * n
		yeildPercent      = 1.0                       // 收率
		geoMeanAccPercent = 1.0                       // 平均准确率

		vSet         = make(map[string]bool)
		vPosRatio    = make(map[string]float64)
		vTypeCounts  = make(map[string]int)
		vTypePercent = make(map[string]float64)
	)
	// 遍历结果
	for sangerPairIndex := range result {
		if Record1Result(primer, result, out, index, sangerPairIndex, offset, vSet) {
			n++
		}
	}
	ln = length * n

	// 按照位置统计
	for key := range vSet {
		spt := strings.Split(key, "-")
		vPosRatio[spt[1]]++
		vTypeCounts[spt[4]]++
	}

	// 计算收率
	for pos := range vPosRatio {
		vPosRatio[pos] /= float64(n)
		yeildPercent *= 1.0 - vPosRatio[pos]
	}

	// 计算几何平均
	geoMeanAccPercent = math.Pow(yeildPercent, 1.0/float64(length)) * 100.0
	yeildPercent *= 100

	// 计算类型比率
	for tp := range vTypeCounts {
		vTypePercent[tp] = float64(vTypeCounts[tp]*100) / float64(ln)
	}

	if n > 0 {
		vErrPercent := float64(len(vSet)*100) / float64(ln)
		vAccPercent := 100.0 - vErrPercent
		resultLine = []any{
			index,
			primer.ID,
			primer.Start + offset, primer.End + offset,
			primer.Start, primer.End,
			length,
			n,
			len(vPosRatio),
			len(vSet),
			vTypeCounts["SNV"],
			vTypeCounts["Insertion"],
			vTypeCounts["Deletion"],
			vTypeCounts["SV"],
			vTypePercent["SNV"],
			vTypePercent["Insertion"],
			vTypePercent["Deletion"],
			vTypePercent["SV"],
			vErrPercent,
			vAccPercent,
			yeildPercent,
			geoMeanAccPercent,
		}
	} else {
		resultLine = []any{
			index,
			primer.ID,
			primer.Start + offset, primer.End + offset,
			primer.Start, primer.End,
			length,
			n,
			len(vPosRatio),
			len(vSet),
			vTypeCounts["SNV"],
			vTypeCounts["Insertion"],
			vTypeCounts["Deletion"],
			vTypeCounts["SV"],
			vTypePercent["SNV"],
			vTypePercent["Insertion"],
			vTypePercent["Deletion"],
			vTypePercent["SV"],
			0.0, 0.0, 0.0, 0.0,
		}
	}
	return
}

func RecordPair(pair *Seq, result map[string][2]*tracy.Result, out *os.File, pairIndex, segOffset int) (resultLines [][]interface{}) {
	var (
		segStart = pair.RefStart + segOffset
		segEnd   = pair.RefEnd + segOffset
	)

	slog.Info("RecordPair", "id", pair.ID, "segStart", segStart, "segEnd", segEnd, "start", pair.Start, "end", pair.End)
	for i, primer := range pair.SubSeq {
		index := fmt.Sprintf("%d.%d", pairIndex, i+1)
		offset := segStart - pair.Start
		resultLine := RecordPrimer(primer, result, out, index, offset)
		resultLines = append(resultLines, resultLine)
	}
	return
}

func RecordSeqPrimer(seq *Seq, result map[string][2]*tracy.Result, prefix string) (resultLines [][]any) {
	out := osUtil.Create(prefix + ".result.txt")
	defer simpleUtil.DeferClose(out)

	for i, pair := range seq.SubSeq {
		lines := RecordPair(pair, result, out, i+1, seq.Start)
		resultLines = append(resultLines, lines...)
	}
	return
}

type PosScore struct {
	Pos int
	Sum float64
	// Key string
}

func GetHightRatio(vPosRatio map[int]float64, k, n int) map[int]float64 {
	// 第一步：构建每个位置的窗口总和
	var scores []PosScore
	for pos := range vPosRatio {
		var sum float64 = 0
		for i := pos - k; i <= pos+k; i++ {
			if val, ok := vPosRatio[i]; ok {
				sum += val
			}
		}
		scores = append(scores, PosScore{Pos: pos, Sum: sum})
	}

	// 第二步：按窗口总和降序排序
	sort.Slice(scores, func(i, j int) bool {
		return scores[i].Sum > scores[j].Sum
	})

	// 第三步：选择非重叠窗口
	used := make(map[int]bool)
	vPosHightRatio := make(map[int]float64)

	for _, score := range scores {
		overlap := false
		for i := score.Pos - k; i <= score.Pos+k; i++ {
			if used[i] {
				overlap = true
				break
			}
		}
		if overlap {
			continue
		}
		if score.Sum >= 4 || score.Sum/float64(n) >= 0.6 {
			// 标记中心位置
			vPosHightRatio[score.Pos] = vPosRatio[score.Pos]
			// 标记这个窗口的所有位置为已使用
			for i := score.Pos - k; i <= score.Pos+k; i++ {
				used[i] = true
			}
		}
	}

	return vPosHightRatio
}

func GetDeletionHightRatio(vPosDeletionRatio map[int]float64, n int) map[int]float64 {
	var vPosDeletionHightRatio = make(map[int]float64)
	for pos := range vPosDeletionRatio {
		ratio := vPosDeletionRatio[pos] / float64(n)
		if ratio >= 0.8 {
			vPosDeletionHightRatio[pos] = ratio
		}
	}
	return vPosDeletionHightRatio
}

func RecordSeq(seq *Seq, result map[string][2]*tracy.Result, prefix string) (resultLine []any) {
	out := osUtil.Create(prefix + ".seq.result.txt")
	defer simpleUtil.DeferClose(out)

	var (
		length            = seq.End - seq.Start // 长度
		n                 = 0                   // 有效Sanger个数
		failN             = 0                   // 无效Sanger个数
		ln                = 0                   // length * n
		yeildPercent      = 1.0                 // 收率
		geoMeanAccPercent = 1.0                 // 平均准确率
		sangerSet         []string

		vSet                   = make(map[string]bool) // 变异
		vPosRatio              = make(map[int]float64) // 变异位置
		vPosHightRatio         = make(map[int]float64) // 高频变异位置
		vPosDeletionRatio      = make(map[int]float64) // 缺失变异位置
		vPosDeletionHightRatio = make(map[int]float64) // 高频缺失
		vTypeCounts            = make(map[string]int)
		vTypePercent           = make(map[string]float64)
		selectClones           = make(map[string]bool) // 正确克隆
		selectHetClones        = make(map[string]bool) // Het正确克隆
	)

	for sangerPairIndex, pairResult := range result {
		if Record1Result(seq, result, out, seq.ID, sangerPairIndex, 0, vSet) {
			n++
			variantCount := 0
			variantHetCount := 0
			sangerSet = append(sangerSet, sangerPairIndex)
			for _, r := range pairResult {
				if r != nil && r.Variants != nil {
					for _, v := range r.Variants.Variants {
						if v.Pos >= seq.Start && v.Pos <= seq.End {
							variantCount++
							if v.Genotype == "het." {
								variantHetCount++
							}
						}
					}
				}
			}
			switch variantCount {
			case 0:
				selectClones[sangerPairIndex] = true
			case variantHetCount: // only het
				selectHetClones[sangerPairIndex] = true
			}
		} else {
			failN++
		}
	}
	ln = length * n

	// 按照位置统计
	// changeID, v.Pos, v.Ref, v.Alt, v.Type
	for key := range vSet {
		spt := strings.Split(key, "-")
		pos := stringsUtil.Atoi(spt[1])
		ref := spt[2]
		alt := spt[3]
		vType := spt[4]
		vPosRatio[pos]++
		vTypeCounts[vType]++
		if spt[4] == "Deletion" {
			for k := range ref {
				if k >= len(alt) {
					vPosDeletionRatio[pos+k]++
				}
			}
		}
	}

	vPosHightRatio = GetHightRatio(vPosRatio, 1, n)
	vPosDeletionHightRatio = GetDeletionHightRatio(vPosDeletionRatio, n)
	// 计算收率
	for pos := range vPosRatio {
		vPosRatio[pos] /= float64(n)
		yeildPercent *= 1.0 - vPosRatio[pos]
	}

	// 计算几何平均
	geoMeanAccPercent = math.Pow(yeildPercent, 1.0/float64(length)) * 100.0
	yeildPercent *= 100

	// 计算类型比率
	for tp := range vTypeCounts {
		vTypePercent[tp] = float64(vTypeCounts[tp]*100) / float64(ln)
	}

	// 正确克隆
	rightCloneIDs := lo.Keys(selectClones)
	rcIDs := strings.Join(rightCloneIDs, "、")
	rcIDs1 := strings.Join(rightCloneIDs[:min(2, len(rightCloneIDs))], "、")
	rcIDs2 := strings.Join(rightCloneIDs[min(2, len(rightCloneIDs)):], "、")
	// Het正确克隆
	rightHetCloneIDs := lo.Keys(selectHetClones)
	rcHetIDs := strings.Join(rightHetCloneIDs, "、")

	vErrPercent := 0.0
	vAccPercent := 0.0
	if n > 0 {
		vErrPercent = float64(len(vSet)*100) / float64(ln)
		vAccPercent = 100.0 - vErrPercent
	} else { // no valid clone
		geoMeanAccPercent = 0.0
		yeildPercent = 0.0
	}
	id1 := seq.ID
	id2 := strings.Replace(id1, "_", "-", 1)
	reg1 := regexp.MustCompile(fmt.Sprintf("^%s-(\\d+)", id1))
	reg2 := regexp.MustCompile(fmt.Sprintf("^%s-(\\d+)", id2))
	if len(rightCloneIDs) > 0 {
		if reg1.MatchString(rightCloneIDs[0]) {
			for _, id := range rightCloneIDs {
				m := reg1.FindStringSubmatch(id)
				if m == nil {
					log.Fatalf("can not parse clone:[%s]vs[%s]", id, seq.ID)
				}
				seq.rIDs = append(seq.rIDs, m[1])
			}
			rcIDs = id1 + "-" + strings.Join(seq.rIDs, "、")
			rcIDs1 = id1 + "-" + strings.Join(seq.rIDs[:min(2, len(seq.rIDs))], "、")
			if len(seq.rIDs) >= 2 {
				rcIDs2 = id1 + "-" + strings.Join(seq.rIDs[min(2, len(seq.rIDs)):], "、")
			}
		} else if reg2.MatchString(rightCloneIDs[0]) {
			for _, id := range rightCloneIDs {
				m := reg2.FindStringSubmatch(id)
				if m == nil {
					log.Fatalf("can not parse clone:[%s]vs[%s]", id, seq.ID)
				}
				seq.rIDs = append(seq.rIDs, m[1])
			}
			rcIDs = id2 + "-" + strings.Join(seq.rIDs, "、")
			rcIDs1 = id2 + "-" + strings.Join(seq.rIDs[:min(2, len(seq.rIDs))], "、")
			if len(seq.rIDs) > 2 {
				rcIDs2 = id2 + "-" + strings.Join(seq.rIDs[min(2, len(seq.rIDs)):], "、")
			}
		} else {
			log.Printf("can not parse clone:[%s]vs[%s|%s]", rightCloneIDs[0], reg1, reg2)
		}
	}

	resultLine = []any{
		seq.ID,
		seq.Start, seq.End,
		length,
		n,
		failN,
		len(selectClones),
		len(selectHetClones),
		len(vPosRatio),
		len(vSet),
		len(vPosHightRatio),
		len(vPosDeletionHightRatio),
		vTypeCounts["SNV"],
		vTypeCounts["Insertion"],
		vTypeCounts["Deletion"],
		vTypeCounts["SV"],
		vTypePercent["SNV"],
		vTypePercent["Insertion"],
		vTypePercent["Deletion"],
		vTypePercent["SV"],
		vErrPercent,
		vAccPercent,
		yeildPercent,
		geoMeanAccPercent,
		sangerSet,
		rcIDs,
		rcIDs1,
		rcHetIDs,
		rcIDs2,
	}
	// 写入 result
	resultFormat := "%s\t%3d-%3d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f%%\t%.4f%%\t%.4f%%\t%.4f%%\t%4.f%%\t%.4f%%\t%.4f%%\t%.4f%%\t%v\t%v\n"
	fmtUtil.Fprintf(out, resultFormat, resultLine...)

	return
}

func CreateFasta(seq *Seq, prefix string) {
	fa := osUtil.Create(prefix + ".fa")
	fmtUtil.Fprintf(fa, ">%s\n%s\n", seq.ID, seq.Seq)
	simpleUtil.CheckErr(fa.Close())
}

func WriteResultTxt(path string, title, lines []string) {
	resultFile := osUtil.Create(path)
	defer simpleUtil.DeferClose(resultFile)

	fmtUtil.FprintStringArray(resultFile, title, "\t")
	fmtUtil.FprintStringArray(resultFile, lines, "\n")
}

func WriteSlice(path, format string, title, list []string, data map[string][][]any) {
	resultFile := osUtil.Create(path)
	defer simpleUtil.DeferClose(resultFile)

	fmtUtil.FprintStringArray(resultFile, title, "\t")
	for _, v := range list {
		for _, s := range data[v] {
			fmtUtil.Fprintf(resultFile, format, s...)
		}
	}
}

func WriteSliceSheet(xlsx *excelize.File, sheet string, title, list []string, data map[string][][]any) {
	simpleUtil.HandleError(xlsx.NewSheet(sheet))
	xlsx.SetSheetRow(sheet, "A1", &title)

	row := 2
	for _, v := range list {
		for _, s := range data[v] {
			xlsx.SetSheetRow(sheet, "A"+strconv.Itoa(row), &s)
			row++
		}
	}
}

func GetTracyStatusLines(id string, result map[string][2]*tracy.Result) (data [][]any) {
	for sangerPairIndex, pairResult := range result {
		for sangerIndex, result := range pairResult {
			slog.Debug("GetTracyStatusLines", "id", id, "sangerPairIndex", sangerPairIndex, "sangerIndex", sangerIndex, "result", result)
			if result == nil {
				continue
			}

			variantCount := 0
			hetCount := 0
			if result.Variants != nil {
				for _, v := range result.Variants.Variants {
					if v.Qual >= MaxQual {
						variantCount++
					}
					if v.Genotype == "het." {
						hetCount++
					}
				}
			}
			row := []any{
				id, sangerPairIndex, sangerIndex,
				result.Status, result.Pass,
				variantCount,
				hetCount,
				result.AlignResult.BoundMatchRatio,
			}
			data = append(data, row)
		}
	}
	return
}

func LoadRawSequence(xlsx *excelize.File, sheet string) (map[string]*Seq, []string) {
	var (
		title []string

		geneList []string
		geneMap  = make(map[string]*Seq)

		rows = simpleUtil.HandleError(xlsx.GetRows(sheet))
	)

	for i, row := range rows {
		if i == 0 {
			title = row
			continue
		}
		var seq = &Seq{}
		for j, cell := range row {
			switch title[j] {
			case "基因名称":
				if cell == "" {
					log.Fatalf("%s [片段名称](%d,%d) 为空", sheet, i+1, j+1)
				}
				seq.ID = cell
			case "DNA序列":
				seq.Seq = cell
			}
		}
		geneMap[seq.ID] = seq
		geneList = append(geneList, seq.ID)
	}
	return geneMap, geneList
}

func LoadSegmentSequence(xlsx *excelize.File, sheet string, geneMap map[string]*Seq) (map[string]*Seq, []string) {
	var (
		title   []string
		seqList []string

		seqMap = make(map[string]*Seq)
		rows   = simpleUtil.HandleError(xlsx.GetRows(sheet))
	)

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
				seq.RefID = cell[:len(cell)-2]
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

		gene, ok := geneMap[seq.RefID]
		if !ok {
			log.Fatalf("can not find Ref[%s] of Segment[%s],[keys:%+v]", seq.RefID, seq.ID, lo.Keys(geneMap))
		}
		gene.SubSeq = append(gene.SubSeq, seq)
	}
	return seqMap, seqList
}

func LoadPrimerPairSequence(xlsx *excelize.File, sheet string, segmentMap map[string]*Seq) {
	var (
		title = []string{}
		rows  = simpleUtil.HandleError(xlsx.GetRows(sheet))
	)
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
		refSeq, ok := segmentMap[pair.RefID]
		if !ok {
			log.Fatalf("ref not found for:[id:%s,ref:%s],[keys:%+v]", pair.ID, pair.RefID, lo.Keys(segmentMap))
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
}

func AddSequencingResultPlate(xlsx *excelize.File, sheet string, geneList []string, geneMap map[string]*Seq) {
	var (
		bgColor1 = "#8CDDFA"
		bgColor2 = "#FFC000"
		bgColor3 = "#FFD966"

		bgStyle1 = &excelize.Style{
			Fill: excelize.Fill{
				Type:    "pattern",
				Color:   []string{bgColor1},
				Pattern: 1,
			},
		}
		bgStyle2 = &excelize.Style{
			Fill: excelize.Fill{
				Type:    "pattern",
				Color:   []string{bgColor2},
				Pattern: 1,
			},
		}
		bgStyle3 = &excelize.Style{
			Fill: excelize.Fill{
				Type:    "pattern",
				Color:   []string{bgColor3},
				Pattern: 1,
			},
		}

		style1 = simpleUtil.HandleError(xlsx.NewStyle(bgStyle1))
		style2 = simpleUtil.HandleError(xlsx.NewStyle(bgStyle2))
		style3 = simpleUtil.HandleError(xlsx.NewStyle(bgStyle3))
	)

	simpleUtil.HandleError(xlsx.NewSheet(sheet))
	xlsx.SetSheetRow(sheet, "A3", &[]string{"基因名称", "长度", "节数"})
	xlsx.SetSheetRow(sheet, "E3", &PlateCols)
	xlsx.SetSheetRow(sheet, "S3", &PlateCols)
	xlsx.SetCellStr(sheet, "D2", "批次号")
	xlsx.SetSheetCol(sheet, "D4", &PlateRows)
	xlsx.SetCellStr(sheet, "R2", "批次号")
	xlsx.SetSheetCol(sheet, "R4", &PlateRows)
	xlsx.SetCellStyle(sheet, "B4", "D11", style1)
	xlsx.SetCellStyle(sheet, "R4", "R11", style1)
	xlsx.SetCellStyle(sheet, "E3", "P3", style2)
	xlsx.SetCellStyle(sheet, "S3", "AD3", style2)

	for i, geneID := range geneList {
		gene := geneMap[geneID]
		xlsx.SetCellStr(sheet, fmt.Sprintf("A%d", 4+i*2), geneID)
		xlsx.SetCellStr(sheet, fmt.Sprintf("B%d", 4+i*2), fmt.Sprintf("%dbp", len(gene.Seq)))
		xlsx.SetCellInt(sheet, fmt.Sprintf("C%d", 4+i*2), len(gene.SubSeq))

		for j, segment := range gene.SubSeq {
			cellName := simpleUtil.HandleError(excelize.CoordinatesToCellName(19+j, 4+i*2))
			xlsx.SetCellInt(sheet, cellName, segment.End-segment.Start)
			cellName = simpleUtil.HandleError(excelize.CoordinatesToCellName(19+j, 4+i*2+1))
			xlsx.SetCellInt(sheet, cellName, segment.End-segment.Start)
			cellName = simpleUtil.HandleError(excelize.CoordinatesToCellName(5+j, 4+i*2))
			xlsx.SetCellStr(sheet, cellName, "N")
			xlsx.SetCellStyle(sheet, cellName, cellName, style3)
			if len(segment.rIDs) > 0 {
				xlsx.SetCellStr(sheet, cellName, segment.ID+"-"+segment.rIDs[0])
			}
			cellName = simpleUtil.HandleError(excelize.CoordinatesToCellName(5+j, 4+i*2+1))
			xlsx.SetCellStr(sheet, cellName, "N")
			xlsx.SetCellStyle(sheet, cellName, cellName, style3)
			if len(segment.rIDs) > 1 {
				xlsx.SetCellStr(sheet, cellName, segment.ID+"-"+strings.Join(segment.rIDs[1:], "、"))
			}
		}
	}
}
