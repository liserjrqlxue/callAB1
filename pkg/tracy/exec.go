package tracy

import (
	"encoding/json"
	"log/slog"
	"os"
	"os/exec"
	"strconv"

	"github.com/liserjrqlxue/goUtil/fmtUtil"
	"github.com/liserjrqlxue/goUtil/osUtil"
	"github.com/liserjrqlxue/goUtil/simpleUtil"
)

func Basecall(tracy, input, prefix string, left, right int) error {
	var args = []string{
		"basecall",
		"-o", prefix + ".basecall.json",
		"-q", strconv.Itoa(left),
		"-u", strconv.Itoa(right),
		input,
	}
	cmd := exec.Command(tracy, args...)
	slog.Info("Basecall", "CMD", cmd)
	cmd.Stderr = os.Stderr
	cmd.Stdout = os.Stdout
	return cmd.Run()
}

func Decompose(tracy, ref, input, prefix string, left, right int) error {
	var args = []string{
		"decompose",
		"-r", ref,
		"-o", prefix + ".decompose",
		"-q", strconv.Itoa(left),
		"-u", strconv.Itoa(right),
		"-v",
		input,
	}
	cmd := exec.Command(tracy, args...)
	slog.Info("Decompose", "CMD", cmd)
	cmd.Stderr = os.Stderr
	cmd.Stdout = os.Stdout
	return cmd.Run()
}

func RunBasecall(tracy, input, prefix string, left, right int) (result Result, err error) {
	err = Basecall(tracy, input, prefix, left, right)
	if err != nil {
		slog.Error("Basecall", "err", err)
		return
	}

	var resultJson []byte
	resultJson, err = os.ReadFile(prefix + ".basecall.json")
	if err != nil {
		slog.Error("Load basecall.json", "json", prefix+".basecall.json", "err", err)
		return
	}

	err = json.Unmarshal(resultJson, &result)
	if err != nil {
		slog.Error("Unmarshal basecall.json", "jsonByte", resultJson, "err", err)
		return
	}

	return
}

func RunDecompose(tracy, ref, input, prefix string, left, right int) (result Result, err error) {
	err = Decompose(tracy, ref, input, prefix, left, right)
	if err != nil {
		slog.Error("Decompose", "err", err)
		return
	}

	var resultJson []byte
	resultJson, err = os.ReadFile(prefix + ".decompose.json")
	if err != nil {
		slog.Error("Load decompose.json", "json", prefix+".decompose.json", "err", err)
		return
	}

	err = json.Unmarshal(resultJson, &result)
	if err != nil {
		slog.Error("Unmarshal decompose.json", "jsonByte", resultJson, "err", err)
		return
	}

	return
}

func RunSingle(tracy, ref, input, prefix string) (result Result, err error) {
	// recover panic
	defer func() {
		if e := recover(); e != nil {
			slog.Error("RunSingle Recover", "e", e, "err", err)
			err = e.(error)
			return
		}
	}()

	var left, right = 0, 0

	result = simpleUtil.HandleError(RunBasecall(tracy, input, prefix, left, right))

	// trim
	left = 60
	right = max(60, len(result.BasecallPos)-700)
	result = simpleUtil.HandleError(RunDecompose(tracy, ref, input, prefix, left, right))

	result.Variants.CalVariants()
	var out = osUtil.Create(prefix + ".decompose.variants.txt")
	defer simpleUtil.DeferClose(out)
	// output title
	fmtUtil.FprintStringArray(out, append(result.Variants.Columns, "xrange1", "xrange2"), "\t")
	fmtUtil.Fprintln(out, result.Variants.String())

	return

}
