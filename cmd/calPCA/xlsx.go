package main

import (
	"github.com/liserjrqlxue/goUtil/simpleUtil"
	"github.com/xuri/excelize/v2"
)

func CoordinatesToCellName(col int, row int, abs ...bool) string {
	return simpleUtil.HandleError(
		excelize.CoordinatesToCellName(
			col, row, abs...,
		),
	)
}
