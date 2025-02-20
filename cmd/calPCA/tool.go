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

func RecordVariant(v *tracy.Variant, out *os.File, index, id string, sangerIndex, start, end, hetCount int, boundMatchRatio float64, variantSet map[string]bool, pass bool) {
	if v.Pos > start && v.Pos <= end {
		fmtUtil.Fprintf(
			out,
			"%s\t%s\t%3d-%3d\t%d.1\t%f\t%d\t%s\n",
			index,
			id,
			start, end,
			sangerIndex,
			boundMatchRatio, hetCount,
			v,
		)
		if pass {
			key := fmt.Sprintf("%d-%d-%s-%s", sangerIndex, v.Pos, v.Ref, v.Alt)
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
		RecordVariant(
			v, out, index, id, sangerIndex, start, end,
			result1.Variants.HetCount,
			result1.AlignResult.BoundMatchRatio,
			variantSet,
			result1.Pass,
		)
	}
	for _, v := range result2.Variants.Variants {
		RecordVariant(
			v, out, index, id, sangerIndex, start, end,
			result2.Variants.HetCount,
			result2.AlignResult.BoundMatchRatio,
			variantSet,
			result2.Pass,
		)
	}
	return
}

func RecordPrimer(primer *Seq, result map[int][2]*tracy.Result, out *os.File, index string, offset int) {
	var (
		length = primer.End - primer.Start // 长度

		n          = 0   // sanger 个数
		yeild      = 1.0 // 收率
		geoMeanAcc = 1.0 // 平均准确率

		variantSet   = make(map[string]bool)
		variantRatio = make(map[string]float64)
	)
	// 遍历结果
	for sangerIndex := range result {
		if Record1Result(primer, result, out, index, sangerIndex, offset, variantSet) {
			n++
		}
	}

	// 按照位置统计
	for key := range variantSet {
		pos := strings.Split(key, "-")[1]
		variantRatio[pos]++
	}

	// 计算收率
	for pos := range variantRatio {
		variantRatio[pos] /= float64(n)
		yeild *= 1.0 - variantRatio[pos]
	}

	// 计算几何平均
	geoMeanAcc = math.Pow(yeild, 1.0/float64(length))

	fmt.Printf(
		"%s\t%s\t%3d-%3d\t%d-%d\t%d\t%d\t%d\t%f\t%f\n",
		index,
		primer.ID,
		primer.Start+offset, primer.End+offset,
		primer.Start, primer.End,
		n,
		len(variantSet),
		len(variantRatio),
		yeild, geoMeanAcc,
	)
}

func RecordPair(pair *Seq, result map[int][2]*tracy.Result, out *os.File, pairIndex, segOffset int) {
	var (
		segStart = pair.RefStart + segOffset
		segEnd   = pair.RefEnd + segOffset
	)

	slog.Info("RecordPair", "id", pair.ID, "segStart", segStart, "segEnd", segEnd, "start", pair.Start, "end", pair.End)
	for i, primer := range pair.SubSeq {
		index := fmt.Sprintf("%d.%d", pairIndex, i+1)
		offset := segStart - pair.Start
		RecordPrimer(primer, result, out, index, offset)
	}
}

func RecordSeq(seq *Seq, result map[int][2]*tracy.Result, prefix string) {
	out := osUtil.Create(prefix + ".result.txt")
	defer simpleUtil.DeferClose(out)

	for i, pair := range seq.SubSeq {
		RecordPair(pair, result, out, i+1, seq.Start)
	}
}

func CreateFasta(seq *Seq, prefix string) {
	fa := osUtil.Create(prefix + ".fa")
	fmtUtil.Fprintf(fa, ">%s\n%s\n", seq.ID, seq.Seq)
	simpleUtil.CheckErr(fa.Close())
}
