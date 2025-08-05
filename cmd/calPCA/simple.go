package main

import (
	"fmt"
	"log"
	"log/slog"
	"os"
	"path/filepath"

	"github.com/liserjrqlxue/DNA/pkg/util"
	"github.com/liserjrqlxue/goUtil/simpleUtil"
	"github.com/xuri/excelize/v2"
)

// 化学补充
func runFix(path string, segmentList []string, segmentMap map[string]*Seq) {
	var (
		xlsx   = excelize.NewFile()
		sheet1 = "化补2清单"
		sheet2 = "化补2引物板位图"

		chemicalSuppList []*util.Seq
	)

	// 化补2清单
	simpleUtil.HandleError(xlsx.NewSheet(sheet1))
	xlsx.SetSheetRow(sheet1, "A1", &[]string{"片段名称", "序列", "长度", "设计成功"})

	for _, id := range segmentList {
		seq, ok := segmentMap[id]
		if ok && seq.CloneHit == 0 {
			id := seq.ID + "C"
			seq := util.NewSeq(id, seq.Seq, filepath.Join(*outputDir, seq.ID, id), true, false)
			chemicalSuppList = append(chemicalSuppList, seq)
		}
	}

	// "化补2引物板位图"
	simpleUtil.HandleError(xlsx.NewSheet(sheet2))
	xlsx.SetSheetCol(sheet2, "A2", &[]string{"金唯智引物订购订单号：", "", "金唯智引物订购单板位图："})
	xlsx.SetSheetCol(sheet2, "B6", &[]string{"A", "B", "C", "D", "E", "F", "G", "H"})
	xlsx.SetSheetRow(sheet2, "C5", &[]int{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12})

	for i, seq := range chemicalSuppList {
		pass := runSimple(seq)

		xlsx.SetSheetRow(sheet1, fmt.Sprintf("A%d", i+2), &[]any{seq.Name, seq.Seq, len(seq.Seq), pass})

		// fmt.Printf("%s\t%s\t%d\n", seq.Name, string(rune('A'+i*2)), seq.PrimerPairCount)
		halfBreak := (seq.PrimerPairCount + 1) / 2
		for j, pair := range seq.SegmentPairs {
			row := string(rune('A' + j/halfBreak*4 + j%halfBreak))
			fmt.Printf(
				"\t%d\t%s%d:%s\t%s%d:%s\n",
				j+1,
				row, i*2+1,
				pair.Left.Name,
				row, i*2+2,
				pair.Right.Name,
			)
			cellName := simpleUtil.HandleError(excelize.CoordinatesToCellName(3+i*2, 6+j/halfBreak*4+j%halfBreak))
			xlsx.SetSheetRow(sheet2, cellName, &[]string{pair.Left.Name, pair.Right.Name})
		}
	}

	simpleUtil.CheckErr(xlsx.DeleteSheet("Sheet1"))
	log.Printf("SaveAs(%s)", path)
	simpleUtil.CheckErr(xlsx.SaveAs(path))
}

func runSimple(seq *util.Seq) (status bool) {
	defer func() {
		if e := recover(); e != nil {
			status = false
		}
	}()
	// slog.SetLogLoggerLevel(slog.LevelDebug)
	// log.Printf("%+v\n", seq)
	// logFile := osUtil.Create("PrimerDesigner.log")
	// slog.SetDefault(slog.New(slog.NewTextHandler(logFile, nil)))
	simpleUtil.CheckErr(os.MkdirAll(filepath.Dir(seq.OutputPrefix), 0755))

	// log
	// util.InitLog(seq)

	// slog.Info("Step1 Done")
	// seq.Step1()
	seq.CalAll()

	if !seq.FindPrimerPairs() || seq.PrimerPairCount > 16 {
		slog.Error("SegmentSplit fail", "seq.Name", seq.Name, "seq.PrimerPairCount", seq.PrimerPairCount)
		return false
	}
	// output after design
	// slog.Info("Step3 Start")
	// seq.Step3()
	return true
}
