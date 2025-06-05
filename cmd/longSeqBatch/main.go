package main

import (
	"flag"
	"log"
	"path/filepath"

	"github.com/liserjrqlxue/version"
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
	version.LogVersion()
	flag.Parse()
	if *inputDir == "" || *result == "" {
		flag.Usage()
		log.Fatal("-i/-o required!")
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
	var genes = &Genes{
		GeneInfo: make(map[string]*Gene),
		GeneList: []string{},
		OutDir:   *outDir,
	}

	genes.LoadInput(*refXlsx)

	genes.GeneRun()

	genes.CreateResult()
}
