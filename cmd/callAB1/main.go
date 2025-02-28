package main

import (
	"callAB1/pkg/tracy"
	"flag"
	"log"
	"path/filepath"
	"strings"
)

// flag
var (
	bin = flag.String(
		"tracy",
		"tracy",
		"https://github.com/gear-genomics/tracy binary path",
	)
	input = flag.String(
		"i",
		"",
		"input ab1 path",
	)
	input2 = flag.String(
		"i2",
		"",
		"paired ab1 path",
	)
	ref = flag.String(
		"ref",
		"",
		"ref fasta path",
	)
	prefix = flag.String(
		"prefix",
		"",
		"output prefix",
	)
	trim = flag.Int(
		"trim",
		60,
		"min trim length",
	)
	maxLength = flag.Int(
		"max",
		700,
		"max length, tail trimed",
	)
)

func init() {
	flag.Parse()
	if *input == "" {
		flag.PrintDefaults()
		log.Fatal("input is required")
	}
	if *ref == "" {
		flag.PrintDefaults()
		log.Fatal("ref is required")
	}
	if *prefix == "" {
		*prefix = strings.TrimSuffix(filepath.Base(*input), filepath.Ext(*input))
	}

	tracy.Trim = *trim
	tracy.MaxLength = *maxLength
}

func main() {

	if *input2 == "" {
		tracy.RunSingle(*bin, *ref, *input, *prefix, false)
	} else {
		tracy.RunPair(*bin, *ref, *input, *input2, *prefix)
	}
}
