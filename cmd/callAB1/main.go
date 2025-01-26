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
		"input",
		"input",
		"input ab1 path",
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
)

func main() {
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

	tracy.RunSingle(*bin, *ref, *input, *prefix)
}
