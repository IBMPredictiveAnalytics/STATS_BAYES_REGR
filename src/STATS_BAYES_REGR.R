#/***********************************************************************
# * Licensed Materials - Property of IBM 
# *
# * IBM SPSS Products: Statistics Common
# *
# * (C) Copyright IBM Corp. 2015
# *
# * US Government Users Restricted Rights - Use, duplication or disclosure
# * restricted by GSA ADP Schedule Contract with IBM Corp. 
# ************************************************************************/

# author__ = "IBM SPSS, JKP"
# version__ = "1.0.0"

# History
# 18-sep-2015 Original Version


gtxt <- function(...) {
    return(gettext(...,domain="STATS_BAYES_REGR"))
}

gtxtf <- function(...) {
    return(gettextf(...,domain="STATS_BAYES_REGR"))
}

kwdmap = list("allmodels"="all", "stepdown"="top", "stepup"="bottom")
### MAIN ROUTINE ###
doBayesregr = function(dep, indep, models="allmodels", comparison=NULL, maxmodels=50000,
    plotbayesf=FALSE, index=NULL, rscalecont="medium", iterations=1000,
    modelsource="none", modelfile=NULL, workspaceaction="clear", modelfileout=NULL) {
    # Estimate Bayes regression
    
    # The modelsource and modelfile
    # parameters are not implemented, awaiting requests for that functionality

    setuplocalization("STATS_BAYES_REGR")
    
    # A warnings proc name is associated with the regular output
    # (and the same omsid), because warnings/errors may appear in
    # a separate procedure block following the regular output
    procname=gtxt("Bayesian Regression")
    warningsprocname = gtxt("Bayesian Regression: Warnings")
    omsid="STATSBAYESREGR"
    warns = Warn(procname=warningsprocname,omsid=omsid)

    tryCatch(library(BayesFactor), error=function(e){
        warns$warn(gtxtf("The R %s package is required but could not be loaded.", "cmprsk"),dostop=TRUE)
        }
    )
    if (!is.null(spssdictionary.GetWeightVariable())) {
        warns$warn(
            gtxt("The dataset is weighted, but case weights are not used in this procedure except for screening out cases with a nonpositive weight"),
            dostop=FALSE)
    }
    if (!is.null(spssdata.GetSplitVariableNames())) {
        warns$warn(
            gtxt("Split variables are not honored by this procedure"),
            dostop=FALSE)
    }
    if (!is.null(comparison) && comparison == 0) {
        comparison = NULL
    }
    # Allow for estimating a single equation
    if (models == "single") {
        comparison = NULL
        index = 1
        plotbayesf = FALSE
    }
    alldata = c(dep, indep)
    frml = paste(dep, paste(indep, collapse="+"), sep="~")
    allargs = as.list(environment())
    dta = spssdata.GetDataFromSPSS(alldata, missingValueToNA=TRUE,
        factorMode="levels")
    if (any(as.logical(lapply(dta,is.factor)))) {
        warns$warn(gtxt("Categorical variables cannot be used in this procedure"),
            dostop=TRUE)
    }
    # The procedure does not allow missing values
    allargs$ncases = nrow(dta)
    dta = dta[complete.cases(dta),]
    allargs$nvalid = nrow(dta)
    if (models != "single") {
        res = tryCatch(regressionBF(as.formula(frml), data=dta, whichModels=kwdmap[models], 
                progress=FALSE, rscaleCont=rscalecont),
            error = function(e) {
                warns$warn(e$message, dostop=TRUE)
            }
        )
    } else {
        res = tryCatch(lmBF(as.formula(frml), data=dta, progress=FALSE,
            rscaleFixed=rscalecont),
            error = function(e) {
                warns$warn(e$message, dostop=TRUE)
            }
        )
    }
    if (!is.null(allargs$comparison)) {
        allargs$comparison = checkcomparison(allargs$comparison, res, warns)

        res = tryCatch(res/res[allargs$comparison],
            error = function(e) {warns$warn(e, dostop=TRUE)}
        )
    }
    
    post = doposterior(allargs, res, warns)
    displayresults(allargs, res, post, warns)
    
    if (!is.null(modelfile)) {
        save(allargsp, res, post, file=modelfile)
    }
    if (workspaceaction == "retain" && is.null(modelfile)) {
        assign("allargsp", allargsest, envir=.GlobalEnv)
        assign("res", res, envir=.GlobalEnv)
        assign("post", post, envir=.GlobalEnv)
    }
    warns$display()
}

checkcomparison = function(comparison, res, warns) {
    # check comparison spec and return if okay
    if (is.null(comparison)) {
        return(NULL)
    }
    if (comparison > length(res)) {
        warns$warn(gtxtf(
            "The comparison or index model number is greater than the number of models, which is %s. Substituting last model", 
            length(res)), dostop=FALSE)
        return(length(res))
    } else {
        return(comparison)
    }
}

doposterior = function(allargs, res, warns) {
    # calculate posterior distribution if model index specified
    
    if (is.null(allargs$index)) {
        return(NULL)
    }
    allargs$index = checkcomparison(allargs$index, res, warns)
    post = posterior(model=res, index=allargs$index, iterations=allargs$iterations,
        progress=FALSE)
    return(post)
}
    
scaletrans=list("medium"=gtxt("medium"), "wide"=gtxt("wide"), "ultrawide"=gtxt("ultrawide"))
waction=list("clear"="clear", "retain"="retain")

displayresults = function(allargs, res, post, warns) {
    # display results
    # allargs is the parameter set
    
    ressum = extractBF(res)

    StartProcedure(allargs[["procname"]], allargs[["omsid"]])
    
    # summary results
    # input specifications
    # although groups can be specified (cengroup), separate results are not
    # produced.
    lbls = c(gtxt("Dependent Variable"),
             gtxt("Comparison Model"),
             gtxt("Number of Cases"),
             gtxt("Number of Valid Cases"),
             gtxt("Prior Scale on Standardized Slopes"),
             gtxt("Posterior Model Index"),
             gtxt("Posterior Iterations"),
             gtxt("Workspace Action"),
             gtxt("Output Model File")
    )

    vals = c(
            allargs$dep,
            ifelse(is.null(allargs$comparison), 
                gtxt("Intercept only"), row.names(ressum)[allargs$comparison]),
            allargs$ncases,
            allargs$nvalid,
            scaletrans[allargs$rscalecont],
            ifelse(is.null(allargs$index), gtxt("--NA--"), allargs$index),
            ifelse(is.null(allargs$index), gtxt("--NA--"), allargs$iterations),
            waction[allargs$workspaceaction],
            ifelse(is.null(allargs$modelfile), gtxt("--NA--"), allargs$modelfile)
    )

    spsspivottable.Display(data.frame(cbind(vals), row.names=lbls), title = gtxt("Summary"),
        collabels=c(gtxt("Summary")), templateName="BAYESREGRSUMMARY", outline=gtxt("Bayes Regression Summary"),
        caption = gtxtf("Computations done by R package BayesFactor, version: %s", packageVersion("BayesFactor"))
    )

    bf = data.frame(seq(1: length(res)),ressum[1:2])
    bf[3] = bf[3] * 100.
    bf = data.frame(bf, length(res) - rank(bf[2]) + 1)
    # add in posterior probabilities excluding Intercept only model
    ###postprob = data.frame(as.BFprobability(newPriorOdds(res) * res))[-(nrow(bf)+1), 1]
    
    # construct posterior probabilities and merge with Bayes factors
    # The order for probabilities may not be the same as for the Bayes factors
    # which requires a few extra steps to get things merged
    # the BF data frame may not have the intercept row, so that row may be discarded
    postprob = data.frame(as.BFprobability(newPriorOdds(res) * res))[1]

    bf = merge(bf, postprob, by="row.names")
    bf = bf[order(bf[[2]]),]
    row.names(bf) = bf[["Row.names"]]
    bf = bf[-1]
    
    ###bf = data.frame(bf, postprob)

    names(bf) = c(gtxt("Model Number"),
        gtxt("Bayes Factor"), gtxt("Error (+-%)"), gtxt("Rank"),
        gtxt("Posterior Probabilities (Equal Prior)"))
    if (allargs$models == "allmodels") {
        caption = gtxt("All models")
    } else if (allargs$models == "stepdown") {
        caption = gtxt("One variable at a time removed from full model")
    } else {
        caption = gtxt("All single variable models")
    }
    spsspivottable.Display(bf,
        title=gtxt("Bayes Factors"),
        rowdim=gtxt("Equation"), 
        hiderowdimtitle=FALSE,
        templateName="BAYESREGRFACTORS",
        outline=gtxt("Bayes Factors"),
        caption=caption
    )
    

    if (!is.null(allargs$index)) {
        postsum = summary(post)
        postsumstats = postsum$statistics[,-4]  # omit time series SEs
        # make sure extras don't conflict with variable names
        row.names(postsumstats)[(nrow(postsumstats)-1): nrow(postsumstats)] = 
            c(gtxt("*Sigma2"), gtxt("*g"))
        names(postsumstats) = c(gtxt("Mean"), gtxt("Std. Deviation"), gtxt("SE Mean"))
        spsspivottable.Display(
            postsumstats, 
            title=gtxtf("Posterior Summary Statistics for Model %s", allargs$index),
            rowdim=gtxt("Variables"),
            hiderowdimtitle=FALSE,
            templateName="BAYESREGRPOSTSTATS",
            outline=gtxt("Posterior Summary Statistics")
        )
        
        postsumquant = postsum$quantiles
        row.names(postsumquant)[(nrow(postsumquant)-1): nrow(postsumquant)] = 
            c(gtxt("*Sigma2"), gtxt("*g"))
        spsspivottable.Display(
            postsumquant,
            title=gtxtf("Posterior Quantiles for Model %s", allargs$index),
            rowdim=gtxt("Variables"),
            hiderowdimtitle=FALSE,
            templateName="BAYESREGRPOSTQUANTILES",
            outline=gtxt("Posterior Quantiles")
        )
    }

    if (allargs$plotbayesf) {
        plot(res)
    }
    
    spsspkg.EndProcedure()
}

Warn = function(procname, omsid) {
    # constructor (sort of) for message management
    lcl = list(
        procname=procname,
        omsid=omsid,
        msglist = list(),  # accumulate messages
        msgnum = 0
    )
    # This line is the key to this approach
    lcl = mylist2env(lcl) # makes this list into an environment

    lcl$warn = function(msg=NULL, dostop=FALSE, inproc=FALSE) {
        # Accumulate messages and, if dostop or no message, display all
        # messages and end procedure state
        # If dostop, issue a stop.

        if (!is.null(msg)) { # accumulate message
            assign("msgnum", lcl$msgnum + 1, envir=lcl)
            # There seems to be no way to update an object, only replace it
            m = lcl$msglist
            m[[lcl$msgnum]] = msg
            assign("msglist", m, envir=lcl)
        } 

        if (is.null(msg) || dostop) {
            lcl$display(inproc)  # display messages and end procedure state
            if (dostop) {
                stop(gtxt("End of procedure"), call.=FALSE)  # may result in dangling error text
            }
        }
    }
    
    lcl$display = function(inproc=FALSE) {
        # display any accumulated messages as a warnings table or as prints
        # and end procedure state, if any

        if (lcl$msgnum == 0) {   # nothing to display
            if (inproc) {
                spsspkg.EndProcedure()
            }
        } else {
            if (!inproc) {
                procok =tryCatch({
                    StartProcedure(lcl$procname, lcl$omsid)
                    TRUE
                    },
                    error = function(e) {
                        FALSE
                    }
                )
            }
            if (procok) {  # build and display a Warnings table if we can
                table = spss.BasePivotTable("Warnings ","Warnings") # do not translate this
                rowdim = BasePivotTable.Append(table,Dimension.Place.row, 
                    gtxt("Message Number"), hideName = FALSE,hideLabels = FALSE)

                for (i in 1:lcl$msgnum) {
                    rowcategory = spss.CellText.String(as.character(i))
                    BasePivotTable.SetCategories(table,rowdim,rowcategory)
                    BasePivotTable.SetCellValue(table,rowcategory, 
                        spss.CellText.String(lcl$msglist[[i]]))
                }
                spsspkg.EndProcedure()   # implies display
            } else { # can't produce a table
                for (i in 1:lcl$msgnum) {
                    print(lcl$msglist[[i]])
                }
            }
        }
    }
    return(lcl)
}

mylist2env = function(alist) {
    env = new.env()
    lnames = names(alist)
    for (i in 1:length(alist)) {
        assign(lnames[[i]],value = alist[[i]], envir=env)
    }
    return(env)
}

# localization initialization
setuplocalization = function(domain) {
    # find and bind translation file names
    # domain is the root name of the extension command .R file, e.g., "SPSSINC_BREUSCH_PAGAN"
    # This would be bound to root location/SPSSINC_BREUSCH_PAGAN/lang

    fpath = Find(file.exists, file.path(.libPaths(), paste(domain, ".R", sep="")))
    bindtextdomain(domain, file.path(dirname(fpath), domain, "lang"))
} 
# override for api to account for extra parameter in V19 and beyond
StartProcedure <- function(procname, omsid) {
    if (substr(spsspkg.GetSPSSVersion(),1, 2) >= 19) {
        spsspkg.StartProcedure(procname, omsid)
    }
    else {
        spsspkg.StartProcedure(omsid)
    }
}



Run = function(args) {
    #Execute the STATS COMPRISK command

    cmdname = args[[1]]
    args = args[[2]]
    oobj = spsspkg.Syntax(list(
        spsspkg.Template("DEPENDENT", subc="", ktype="existingvarlist", var="dep"),
        spsspkg.Template("INDEP", subc="", ktype="existingvarlist", var="indep", islist=TRUE),
        spsspkg.Template("MODELS", subc="", ktype="str", var="models",
            vallist=list("allmodels", "stepdown", "stepup", "single")),
        spsspkg.Template("MODELSOURCE", subc="", ktype="str", var="modelsource",
            vallist=list("none", "workspace", "modelfile")),
        spsspkg.Template("MODELFILE", subc="", ktype="literal", var="modelfile"),
        
        spsspkg.Template("COMPARISON", subc="OPTIONS", ktype="int", var="comparison"),
        spsspkg.Template("MAXMODELS", subc="OPTIONS", ktype="str", var="maxmodels"),
        spsspkg.Template("PLOTMODELS", subc="OPTIONS", ktype="bool", var="plotbayesf"),
        spsspkg.Template("POSTERIORINDEX", subc="OPTIONS", ktype="int", var="index"),
        spsspkg.Template("ITERATIONS", subc="OPTIONS", ktype="int", var="iterations",
            vallist=list(2)),
        spsspkg.Template("PRIORSCALE", subc="OPTIONS", ktype="str", var="rscalecont",
            vallist=list("medium", "wide", "ultrawide")),
        
        spsspkg.Template("WORKSPACE", subc="SAVE", ktype="str", var="workspaceaction",
            vallist=list("retain", "clear")),
        spsspkg.Template("MODELFILE", subc="SAVE", ktype="literal", var="modelfileout")
    ))

    # A HELP subcommand overrides all else
    if ("HELP" %in% attr(args,"names")) {
        helper(cmdname)
    }
    else {
        res <- spsspkg.processcmd(oobj, args, "doBayesregr")
    }
}

helper = function(cmdname) {
    # find the html help file and display in the default browser
    # cmdname may have blanks that need to be converted to _ to match the file
    
    fn = gsub(" ", "_", cmdname, fixed=TRUE)
    thefile = Find(file.exists, file.path(.libPaths(), fn, "markdown.html"))
    if (is.null(thefile)) {
        print("Help file not found")
    } else {
        browseURL(paste("file://", thefile, sep=""))
    }
}
if (exists("spsspkg.helper")) {
assign("helper", spsspkg.helper)
}
