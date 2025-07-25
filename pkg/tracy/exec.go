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

func Basecall(tracy, input, prefix string, left, right int, stdout, stderr *os.File) error {
	var args = []string{
		"basecall",
		"-o", prefix + ".basecall.json",
		"-q", strconv.Itoa(left),
		"-u", strconv.Itoa(right),
		input,
	}
	cmd := exec.Command(tracy, args...)
	slog.Info("Basecall", "CMD", cmd)
	cmd.Stderr = stderr
	cmd.Stdout = stdout
	return cmd.Run()
}

func Align(tracy, ref, input, prefix string, left, right int, stdout, stderr *os.File) error {
	var args = []string{
		"align",
		"-r", ref,
		"-o", prefix + ".align",
		"-q", strconv.Itoa(left),
		"-u", strconv.Itoa(right),
		input,
	}
	cmd := exec.Command(tracy, args...)
	slog.Info("Basecall", "CMD", cmd)
	cmd.Stderr = stderr
	cmd.Stdout = stdout
	return cmd.Run()
}

func Decompose(tracy, ref, input, prefix string, left, right int, stdout, stderr *os.File) error {
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
	cmd.Stderr = stderr
	cmd.Stdout = stdout
	return cmd.Run()
}

func RunBasecall(tracy, input, prefix string, left, right int, stdout, stderr *os.File, override bool) (result Result, err error) {
	jsonFile := prefix + ".basecall.json"
	if override || !osUtil.FileExists(jsonFile) {
		err = Basecall(tracy, input, prefix, left, right, stdout, stderr)
		if err != nil {
			slog.Error("Basecall", "err", err)
			return
		}
	}

	var resultJson []byte
	resultJson, err = os.ReadFile(jsonFile)
	if err != nil {
		slog.Error("Load basecall.json", "json", jsonFile, "err", err)
		return
	}

	err = json.Unmarshal(resultJson, &result)
	if err != nil {
		slog.Error("Unmarshal basecall.json", "jsonByte", resultJson, "err", err)
		return
	}

	return
}

func RunAlign(tracy, ref, input, prefix string, left, right int, stdout, stderr *os.File, override bool) (result AlignResult, err error) {
	jsonFile := prefix + ".align.json"
	if override || !osUtil.FileExists(jsonFile) {
		err = Align(tracy, ref, input, prefix, left, right, stdout, stderr)
		if err != nil {
			slog.Error("Align", "err", err)
			return
		}
	}

	var resultJson []byte
	resultJson, err = os.ReadFile(prefix + ".align.json")
	if err != nil {
		slog.Error("Load align.json", "json", prefix+".align.json", "err", err)
		return
	}

	err = json.Unmarshal(resultJson, &result)
	if err != nil {
		slog.Error("Unmarshal align.json", "jsonByte", resultJson, "err", err)
		return
	}

	return
}

func RunDecompose(tracy, ref, input, prefix string, left, right int, stdout, stderr *os.File, override bool) (result Result, err error) {
	jsonFile := prefix + ".decompose.json"
	if override || !osUtil.FileExists(jsonFile) {
		err = Decompose(tracy, ref, input, prefix, left, right, stdout, stderr)
		if err != nil {
			slog.Error("Decompose", "err", err)
			result = Result{Status: "DecomposeFail", Pass: false}
			os.WriteFile(jsonFile, simpleUtil.HandleError(json.Marshal(result)), 0644)
			return
		}
	}

	var resultJson []byte
	resultJson, err = os.ReadFile(prefix + ".decompose.json")
	if err != nil {
		slog.Error("Load decompose.json", "json", prefix+".decompose.json", "err", err)
		return Result{Status: "DecomposeLoadFail", Pass: false}, err
	}

	err = json.Unmarshal(resultJson, &result)
	if err != nil {
		slog.Error("Unmarshal decompose.json", "jsonByte", resultJson, "err", err)
		return Result{Status: "DecomposeUnmarshalFail", Pass: false}, err
	}

	return
}

func RunSingle(tracy, ref, input, prefix string, override bool) (result Result, err error) {
	// recover panic
	defer func() {
		if e := recover(); e != nil {
			slog.Error("RunSingle Recover", "e", e, "err", err)
			err = e.(error)
			return
		}
	}()

	var left, right = 0, 0

	var stdout = osUtil.Create(prefix + ".stdout.txt")
	defer simpleUtil.DeferClose(stdout)
	var stderr = osUtil.Create(prefix + ".stderr.txt")
	defer simpleUtil.DeferClose(stderr)

	result = simpleUtil.HandleError(RunBasecall(tracy, input, prefix, left, right, stdout, stderr, override))

	// trim
	// left = Trim
	// right = max(Trim, len(result.BasecallPos)-MaxLength)

	alignResult := simpleUtil.HandleError(RunAlign(tracy, ref, input, prefix, left, right, stdout, stderr, override))
	alignResult.CalAlign()
	var summary = osUtil.Create(prefix + ".align.summary.txt")
	defer simpleUtil.DeferClose(summary)
	alignResult.Summary(summary)

	result, err = RunDecompose(tracy, ref, input, prefix, left, right, stdout, stderr, override)
	result.AlignResult = &alignResult
	if err == nil && result.Variants != nil {
		result.Variants.CalVariants()
		var out = osUtil.Create(prefix + ".decompose.variants.txt")
		defer simpleUtil.DeferClose(out)
		// output title
		fmtUtil.FprintStringArray(out, append(result.Variants.Columns, "xrange1", "xrange2"), "\t")
		fmtUtil.Fprintln(out, result.Variants.String())
	} else {
		slog.Error("RunDecompose", "err", err)
	}

	// 评估
	if result.AlignResult.BoundMatchRatio < BoundMatchRaitoLimit {
		result.Status += "LowMatch"
	}

	// if result.Variants != nil && result.Variants.HetCount > HetCountLimit {
	// 	result.Status += "HighHet"
	// }
	if len(result.BasecallPos) < 100 {
		result.Status += "TooShort"
	}
	if result.Status == "" {
		result.Pass = true
		result.Status = "PASS"
	}
	slog.Info("RunSingle", "Status", result.Status)

	return
}

func RunPair(tracy, ref, input1, input2, prefix string) (result1, result2 Result, err error) {
	// recover panic
	defer func() {
		if e := recover(); e != nil {
			slog.Error("RunPair Recover", "e", e, "err", err)
			err = e.(error)
			return
		}
	}()

	result1, err = RunSingle(tracy, ref, input1, prefix+"_1", false)
	if err != nil {
		slog.Error("RunSingle 1", "err", err)
	}
	result2, err = RunSingle(tracy, ref, input2, prefix+"_2", false)
	if err != nil {
		slog.Error("RunSingle 2", "err", err)
	}

	var out = osUtil.Create(prefix + ".decompose.variants.txt")
	defer simpleUtil.DeferClose(out)
	// output title
	var title = []string{"Source"}
	title = append(title, result1.Variants.Columns...)
	title = append(title, "xrange1", "xrange2")
	fmtUtil.FprintStringArray(out, title, "\t")
	for _, v := range result1.Variants.Variants {
		fmtUtil.Fprintf(out, "%s\t%s\n", "1", v.String())
	}
	for _, v := range result2.Variants.Variants {
		fmtUtil.Fprintf(out, "%s\t%s\n", "2", v.String())
	}

	return

}
