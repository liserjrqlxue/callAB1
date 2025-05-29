package main

import (
	"flag"
	"log"
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

		genes = &Genes{
			GeneInfo: make(map[string]*Gene),
			GeneList: []string{},
			OutDir:   *outDir,
		}
	)

	sheet1Data := GetRows2MapArray(xlsx, "Sheet1")
	genes.GetGenes(sheet1Data)

	sheet2Data := GetRows2MapArray(xlsx, "Sheet2")
	genes.LoadSangerFromGlob(sheet2Data)

	genes.GeneRun()

	genes.CreateResult()
}
