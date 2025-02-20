package main

import (
	"callAB1/pkg/tracy"
	"fmt"
	"log"
	"log/slog"
	"path/filepath"
	"strings"

	"github.com/liserjrqlxue/goUtil/osUtil"
	"github.com/liserjrqlxue/goUtil/simpleUtil"
)

func RunTracyBatch(id, prefix, bin string) map[int][2]*tracy.Result {
	// 遍历分析sanger文件 -> result
	result := make(map[int][2]*tracy.Result)
	tag := []byte(strings.Split(id, "_")[1][:4])
	tag[2] = 'B'
	for j := 1; j < 33; j++ {
		prefix_j := fmt.Sprintf("%s_%d", prefix, j)
		path1 := filepath.Join(*sangerDir, fmt.Sprintf("BGA-%s-TY1-%d-T7.ab1", tag, j))
		path2 := filepath.Join(*sangerDir, fmt.Sprintf("BGA-%s-TY1-%d-T7-Term.ab1", tag, j))

		ok1 := osUtil.FileExists(path1)
		ok2 := osUtil.FileExists(path2)
		if ok1 && ok2 {
			// run tracy
			result1, result2, err := tracy.RunPair(bin, prefix+".fa", path1, path2, prefix_j)
			simpleUtil.CheckErr(err)
			result[j] = [2]*tracy.Result{&result1, &result2}
		} else {
			slog.Debug("sanger no found", "id", id, "j", j, "path1", path1, "ok1", ok1, "path2", path2, "ok2", ok2)
			if !ok1 && !ok2 {
				continue
			} else {
				log.Fatalf("sanger pair fail:[%s:%v,%s:%v]", path1, ok1, path2, ok2)
			}
		}
	}

	return result
}
