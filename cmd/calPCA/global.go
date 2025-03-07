package main

var (
	SeqTitle = []string{
		"片段名称",
		"片段范围-起始",
		"片段范围-终止",
		"有效长度",
		"有效Sanger个数",
		"无效Sanger个数",
		"正确克隆个数",
		"变异位置数",
		"变异个数",
		"SNV",
		"INS",
		"DEL",
		"SV",
		"SNV比率(%)",
		"INS比率(%)",
		"DEL比率(%)",
		"SV比率(%)",
		"变异比率(%)",
		"准确率(%)",
		"参考收率(%)",
		"参考单步准确率(%)",
	}

	ResultTitle = []string{
		"index",
		"引物名称",
		"片段范围-起始",
		"片段范围-终止",
		"引物对范围-起始",
		"引物对范围-终止",
		"引物有效长度",
		"有效Sanger个数",
		"变异位置数",
		"变异个数",
		"SNV",
		"INS",
		"DEL",
		"SV",
		"SNV比率(%)",
		"INS比率(%)",
		"DEL比率(%)",
		"SV比率(%)",
		"变异比率(%)",
		"准确率(%)",
		"参考收率(%)",
		"参考单步准确率(%)",
	}
	CloneVariantTitle = []string{
		"Chr",
		"Pos",
		"Ref",
		"Alt",
		"Type",
		"CloneID",
		"SangerCount",
	}
	SetVariantTitle = []string{
		"Chr",
		"Pos",
		"Ref",
		"Alt",
		"Type",
		"CloneCount",
		"ClonePass",
		"CloneRatio",
	}

	tracyStatusTitle = []string{
		"序列名称",
		"克隆号",
		"Sanger",
		"状态",
		"PASS",
		"变异个数",
		"杂合变异个数",
		"比对率",
	}

	CloneCountLimit = 32
)
