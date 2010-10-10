.onLoad <- function(libname, pkgname) {
	library(RGtk2)
		
	if (!is.null(try(gtkCheckVersion(2,16,0),TRUE))) {
		warning("This package needs gtk+ with version 2.16 or higher. ")
	}
	gtkListStoreGetType()
	gtkWindowGetType()
	gtkNotebookGetType()
	gtkVBoxGetType()
	gtkDrawingAreaGetType()
	gtkHBoxGetType()
	gtkLabelGetType()
	gtkEntryGetType()
	gtkHButtonBoxGetType()
	gtkButtonGetType()
	gtkHScaleGetType()
	gtkSpinButtonGetType()
	gtkHButtonBoxGetType()
	gtkVButtonBoxGetType()
	gtkTableGetType()
	gtkCheckButtonGetType()
	gtkScrolledWindowGetType()
	gtkRadioButtonGetType()
	gtkTreeViewGetType()
	gtkTreeViewColumnGetType()
	gtkCellRendererTextGetType()
	gtkAdjustmentGetType()
	gtkViewportGetType()

}
