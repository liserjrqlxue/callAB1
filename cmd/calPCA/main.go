package main

import (
	"callAB1/pkg/tracy"
	"flag"
	"fmt"
	"log"
	"math"
	"os"
	"path/filepath"
	"strings"

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

	for i, id := range seqList {
		seq := seqMap[id]
		fmt.Printf("%d\t%s\t\t%3d-%-3d\n", i+1, seq.ID, seq.Start, seq.End)

		// 创建 Fasta
		prefix := filepath.Join(*outputDir, id)
		fa := osUtil.Create(prefix + ".fa")
		fmtUtil.Fprintf(fa, ">%s\n%s\n", id, seq.Seq)
		simpleUtil.CheckErr(fa.Close())

		resultF := osUtil.Create(prefix + ".result.txt")

		tag := []byte(strings.Split(id, "_")[1][:4])
		tag[2] = 'B'
		result := make(map[int][2]*tracy.Result)
		for j := 1; j < 33; j++ {
			prefix_j := fmt.Sprintf("%s_%d", prefix, j)
			path1 := filepath.Join(*sangerDir, fmt.Sprintf("BGA-%s-TY1-%d-T7.ab1", tag, j))
			path2 := filepath.Join(*sangerDir, fmt.Sprintf("BGA-%s-TY1-%d-T7-Term.ab1", tag, j))

			ok1 := osUtil.FileExists(path1)
			ok2 := osUtil.FileExists(path2)
			if ok1 && ok2 {
				// run tracy
				result1, result2, err := tracy.RunPair(*bin, prefix+".fa", path1, path2, prefix_j)
				simpleUtil.CheckErr(err)
				result[j] = [2]*tracy.Result{&result1, &result2}
			} else if !ok1 && !ok2 {
				log.Printf("sanger no found:[%s:%v,%s:%v]", path1, ok1, path2, ok2)
				continue
			} else {
				log.Fatalf("sanger pair fail:[%s:%v,%s:%v]", path1, ok1, path2, ok2)
			}
		}

		for j, pair := range seq.SubSeq {
			fmt.Printf("%d.%d\t%s\t\t%3d-%3d\t%d\t%d\n", i+1, j+1, pair.ID, pair.RefStart+seq.Start, pair.RefEnd+seq.Start, pair.Start, pair.End)
			for k, primer := range pair.SubSeq {
				start := primer.Start - pair.Start + pair.RefStart + seq.Start
				end := primer.End - pair.Start + pair.RefStart + seq.Start
				length := end - start

				variantSet := make(map[string]bool)
				fmt.Printf(
					"%d.%d.%d\t%s\t%3d-%3d\t%d\t%d\n",
					i+1, j+1, k+1,
					primer.ID,
					start, end,
					primer.Start, primer.End,
				)
				n := 0
				// 遍历结果
				for l := range result {
					keep := false

					boundMatchRatio1 := result[l][0].AlignResult.BoundMatchRatio
					hetCount1 := result[l][0].Variants.HetCount
					boundMatchRatio2 := result[l][1].AlignResult.BoundMatchRatio
					hetCount2 := result[l][1].Variants.HetCount
					for _, v := range result[l][0].Variants.Variants {
						if v.Pos > start && v.Pos <= end {
							fmtUtil.Fprintf(
								resultF,
								"%d.%d.%d\t%s\t%3d-%3d\t%d.1\t%f\t%d\t%s\n",
								i+1, j+1, k+1,
								primer.ID,
								start, end,
								l,
								boundMatchRatio1, hetCount1,
								v,
							)
							if boundMatchRatio1 > 0.9 && hetCount1 < 50 {
								keep = true
								key := fmt.Sprintf("%d-%d-%s-%s", l, v.Pos, v.Ref, v.Alt)
								variantSet[key] = true
							}
						}
					}
					for _, v := range result[l][1].Variants.Variants {
						if v.Pos > start && v.Pos <= end {
							fmtUtil.Fprintf(
								resultF,
								"%d.%d.%d\t%s\t%3d-%3d\t%d.2\t%f\t%d\t%s\n",
								i+1, j+1, k+1,
								primer.ID,
								start, end,
								l,
								boundMatchRatio2, hetCount2,
								v,
							)
							if boundMatchRatio2 > 0.9 && hetCount2 < 50 {
								keep = true
								key := fmt.Sprintf("%d-%d-%s-%s", l, v.Pos, v.Ref, v.Alt)
								variantSet[key] = true
							}
						}
					}
					if keep {
						n++
					}
				}
				// 按照位置统计
				variantRatio := make(map[string]float64)
				for key := range variantSet {
					pos := strings.Split(key, "-")[1]
					variantRatio[pos]++
				}
				// 计算收率
				yeild := 1.0
				for pos := range variantRatio {
					variantRatio[pos] /= float64(n)
					yeild *= 1.0 - variantRatio[pos]
				}
				// 计算几何平均
				geoMeanAcc := math.Pow(yeild, 1.0/float64(length))
				fmt.Printf("%d.%d.%d\t%s\t%3d-%3d\t%d\t%f\n", i+1, j+1, k+1, primer.ID, start, end, length, geoMeanAcc)
			}
		}

		simpleUtil.CheckErr(resultF.Close())

	}
}
