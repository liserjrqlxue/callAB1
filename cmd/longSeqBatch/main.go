package main

import (
	"flag"
	"log"
	"os"
	"path/filepath"

	"github.com/liserjrqlxue/goUtil/simpleUtil"
	"github.com/xuri/excelize/v2"
)

// flag

var (
	bin = flag.String(
		"tracy",
		"tracy",
		"https://github.com/gear-genomics/tracy binary path",
	)
	inputDir = flag.String(
		"i",
		"",
		"input dir",
	)
	refXlsx = flag.String(
		"r",
		"目标序列.xlsx",
		"输入信息",
	)
	seqDir = flag.String(
		"s",
		"测序目录",
		"ab1 dir",
	)
	outDir = flag.String(
		"d",
		"",
		"output dir, default same as -i",
	)
	result = flag.String(
		"o",
		"result.xlsx",
		"result file name",
	)
)

func init() {
	flag.Parse()
	if *inputDir == "" || *result == "" {
		flag.Usage()
		log.Fatal("-i/-o required")
	}
	if *outDir == "" {
		*outDir = *inputDir
	}

	*refXlsx = filepath.Join(*inputDir, *refXlsx)
	*seqDir = filepath.Join(*inputDir, *seqDir)
	*result = filepath.Join(*outDir, *result)
}

var (
	TermTrim = 20
)

func main() {
	var (
		xlsx = simpleUtil.HandleError(excelize.OpenFile(*refXlsx))

		geneInfo map[string]*Gene
		geneList []string
	)
	sheet1Data := GetRows2MapArray(xlsx, "Sheet1")
	geneInfo, geneList = CreateGeneInfoFromDataArray(sheet1Data, *outDir)

	sheet2Data := GetRows2MapArray(xlsx, "Sheet2")
	LoadSangerFromDataArray(geneInfo, sheet2Data)

	simpleUtil.CheckErr(os.MkdirAll(filepath.Join(*outDir, "ref"), 0755))

	for _, geneID := range geneList {
		log.Printf("loop GeneID:[%s]", geneID)
		gene := geneInfo[geneID]
		simpleUtil.CheckErr(gene.CreateRef())

		for cloneID := range gene.Clones {
			CloneRun(gene, cloneID)
		}
	}

	CreateResult(geneInfo, geneList)
}
