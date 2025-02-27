package main

import (
	"callAB1/pkg/tracy"
	"fmt"
	"log"
	"log/slog"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"strings"

	"github.com/liserjrqlxue/goUtil/fmtUtil"
	"github.com/liserjrqlxue/goUtil/osUtil"
	"github.com/liserjrqlxue/goUtil/simpleUtil"
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
	result1, err := tracy.RunSingle(bin, prefix+".fa", path, sangerPrefix)
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

	for _, path := range files {
		RunTracyCY0130(id, prefix, bin, path, result)
	}

	return result
}

func RecordVariant(v *tracy.Variant, out *os.File, index, id, sangerPairIndex string, sangerIndex, start, end, hetCount int, boundMatchRatio float64, variantSet map[string]bool, pass bool) {
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
			key := fmt.Sprintf("%s-%d-%s-%s-%s", sangerPairIndex, v.Pos, v.Ref, v.Alt, v.Type)
			variantSet[key] = true
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

func RecordPrimer(primer *Seq, result map[string][2]*tracy.Result, out *os.File, index string, offset int) (resultLine []interface{}) {
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
		resultLine = []interface{}{
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

func RecordSeq(seq *Seq, result map[string][2]*tracy.Result, prefix string) (resultLines [][]interface{}) {
	out := osUtil.Create(prefix + ".result.txt")
	defer simpleUtil.DeferClose(out)

	for i, pair := range seq.SubSeq {
		lines := RecordPair(pair, result, out, i+1, seq.Start)
		resultLines = append(resultLines, lines...)
	}
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

func WriteSlice(path, format string, title, list []string, data map[string][][]interface{}) {
	resultFile := osUtil.Create(path)
	defer simpleUtil.DeferClose(resultFile)

	fmtUtil.FprintStringArray(resultFile, title, "\t")
	for _, v := range list {
		for _, s := range data[v] {
			fmtUtil.Fprintf(resultFile, format, s...)
		}
	}
}

func GetTracyStatusLines(id string, result map[string][2]*tracy.Result) (data [][]interface{}) {
	for sangerPairIndex, pairResult := range result {
		for sangerIndex, result := range pairResult {
			slog.Debug("GetTracyStatusLines", "id", id, "sangerPairIndex", sangerPairIndex, "sangerIndex", sangerIndex, "result", result)
			if result == nil {
				continue
			}

			variantCount := 0
			hetCount := 0
			if result.Variants != nil {
				variantCount = len(result.Variants.Variants)
				hetCount = result.Variants.HetCount
			}
			row := []interface{}{
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
