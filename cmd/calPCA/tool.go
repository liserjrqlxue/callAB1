package main

import (
	"callAB1/pkg/tracy"
	"fmt"
	"log"
	"log/slog"
	"math"
	"os"
	"path/filepath"
	"strings"

	"github.com/liserjrqlxue/goUtil/fmtUtil"
	"github.com/liserjrqlxue/goUtil/osUtil"
	"github.com/liserjrqlxue/goUtil/simpleUtil"
)

func RunTracy(id, tag, prefix, bin string, sangerIndex int, result map[int][2]*tracy.Result) {
	sangerPrefix := fmt.Sprintf("%s_%d", prefix, sangerIndex)
	path1 := filepath.Join(*sangerDir, fmt.Sprintf("BGA-%s-TY1-%d-T7.ab1", tag, sangerIndex))
	path2 := filepath.Join(*sangerDir, fmt.Sprintf("BGA-%s-TY1-%d-T7-Term.ab1", tag, sangerIndex))

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

// 遍历分析sanger文件 -> result
func RunTracyBatch(id, prefix, bin string) map[int][2]*tracy.Result {
	result := make(map[int][2]*tracy.Result)
	tag := []byte(strings.Split(id, "_")[1][:4])
	tag[2] = 'B'
	for sangerIndex := 1; sangerIndex < 33; sangerIndex++ {
		RunTracy(id, string(tag), prefix, bin, sangerIndex, result)
	}
	return result
}

func RecordVariant(v *tracy.Variant, out *os.File, index, id string, sangerPairIndex, sangerIndex, start, end, hetCount int, boundMatchRatio float64, variantSet map[string]bool, pass bool) {
	if v.Pos > start && v.Pos <= end {
		fmtUtil.Fprintf(
			out,
			"%s\t%s\t%3d-%3d\t%d.%d\t%f\t%d\t%v\t%s\n",
			index,
			id,
			start, end,
			sangerPairIndex, sangerIndex,
			boundMatchRatio, hetCount,
			pass,
			v,
		)
		if pass {
			key := fmt.Sprintf("%d-%d-%s-%s-%s", sangerPairIndex, v.Pos, v.Ref, v.Alt, v.Type)
			variantSet[key] = true
		}
	}
}

func Record1Result(primer *Seq, result map[int][2]*tracy.Result, out *os.File, index string, sangerIndex, offset int, variantSet map[string]bool) (keep bool) {
	var (
		id    = primer.ID
		start = primer.Start + offset
		end   = primer.End + offset

		result1 = result[sangerIndex][0]
		result2 = result[sangerIndex][1]
	)

	if result1.Pass || result2.Pass {
		keep = true
	}

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
	for _, v := range result2.Variants.Variants {
		RecordVariant(
			v, out, index, id, sangerIndex, 2, start, end,
			result2.Variants.HetCount,
			result2.AlignResult.BoundMatchRatio,
			variantSet,
			result2.Pass,
		)
	}
	return
}

func RecordPrimer(primer *Seq, result map[int][2]*tracy.Result, out *os.File, index string, offset int) (resultLine []interface{}) {
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

func RecordPair(pair *Seq, result map[int][2]*tracy.Result, out *os.File, pairIndex, segOffset int) (resultLines [][]interface{}) {
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

func RecordSeq(seq *Seq, result map[int][2]*tracy.Result, prefix string) (resultLines [][]interface{}) {
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

func WriteSlice(path, format string, title []string, data [][]interface{}) {
	resultFile := osUtil.Create(path)
	defer simpleUtil.DeferClose(resultFile)

	fmtUtil.FprintStringArray(resultFile, title, "\t")
	for _, v := range data {
		fmtUtil.Fprintf(resultFile, format, v...)
	}
}

func GetTracyStatusLines(id string, result map[int][2]*tracy.Result) (data [][]interface{}) {
	for sangerPairIndex, pairResult := range result {
		for sangerIndex, result := range pairResult {
			row := []interface{}{
				id, sangerPairIndex, sangerIndex,
				result.Status, result.Pass,
				len(result.Variants.Variants),
				result.Variants.HetCount,
				result.AlignResult.BoundMatchRatio,
			}
			data = append(data, row)
		}
	}
	return
}
