.uiFile <- function(name,package) {
	if (Sys.getenv("LANGUAGE") == "") {
		if (length(grep("French",Sys.getlocale("LC_CTYPE")))>0) {
			name <- paste(name,"_FR",sep="")
		} else {
			name <- paste(name,"_EN",sep="")
		}	
	} else if (Sys.getenv("LANGUAGE") == "fr") {
		name <- paste(name,"_FR",sep="")
	} else {
		name <- paste(name,"_EN",sep="")
	}
	name <- paste(name,".glade",sep="")
	uiFile <- system.file("ui",name,package=package)
	uiFile
}

