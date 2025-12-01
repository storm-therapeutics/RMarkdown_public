library(ggplot2)
library(tidyr)# for 'pivot_wider'
library(magrittr) # for '%<>%'
library(clusterProfiler) # for 'dotplot'
library(pheatmap)
library(readxl) # for 'read_excel'


#' Prepare sample information table
#'
#' @param inventory.file Path to the inventory file (CSV)
#' @param use.time Use time (not just date) information when calculating time differences?
prepare.sample.info <- function(inventory.file, use.time=FALSE) {
  sample.info <- read.csv(inventory.file, check.names=FALSE)
  ## remove unnecessary rows/columns:
  if (any(duplicated(sample.info[["Accession ID"]]))) {
    ind <- sample.info[["Specimen Type Received"]] == "Peripheral Blood-Elution Tube"
    sample.info <- sample.info[ind, ]
    sample.info <- sample.info[!duplicated(sample.info[["Accession ID"]]), ]
  }
  sample.info <- sample.info[, c("Accession ID", "Subject ID", "Visit ID", "Collection Date",
                                 "Collection Time", "Year (DOB)", "Nucleic Acid Concentration (MOL)")]
  ## order by patient and date:
  sample.info[["Collection Date"]] %<>% as.Date(format="%d-%b-%y")
  ord <- order(sample.info[["Subject ID"]], sample.info[["Collection Date"]], sample.info[["Collection Time"]])
  sample.info <- sample.info[ord, ]
  ## add day number (or full time difference):
  parts <- split(sample.info, sample.info[["Subject ID"]])
  parts <- lapply(parts, function(part) {
    ind <- which(part[["Visit ID"]] == "Cycle 1 Day 1")
    ## date of initial visit may be missing for patients that were split across batches:
    if (length(ind) == 0) {
      part[["Time Point"]] <- NA
    }
    else {
      stopifnot(length(ind) == 1)
      if (use.time) {
        times <- part[["Collection Time"]]
        ## impute missing times ("12:00" is a better guess than "0:00"):
        times[times == ""] <- "12:00"
        times <- as.POSIXct(paste(part[["Collection Date"]], times))
        time0 <- times[ind]
        part[["Time Point"]] <- times - time0
        units(part[["Time Point"]]) <- "days"
      } else {
        day0 <- part[ind, "Collection Date"]
        part[["Time Point"]] <- part[["Collection Date"]] - day0
      }
    }
    if (any(part[["Time Point"]] < 0))
      warning("Negative time difference observed for patient ", part[1, "Subject ID"])
    part
  })
  do.call(rbind, parts)
}


#' Match sample accessions to RCC files
#'
#' @param accessions Vector of accessions
#' @param rcc.files Vector of RCC file paths
#' @param allow.unmatched Is it an error if not all accessions can be matched?
match.rcc.files <- function(accessions, rcc.files, allow.unmatched=FALSE) {
  matches <- sapply(accessions, function(id) {
    ind <- grep(id, rcc.files, fixed=TRUE)
    if (length(ind) == 0) {
      if (allow.unmatched) return(NA)
      stop("no RCC file for accession ", id)
    }
    if (length(ind) > 1) stop("multiple RCC files for accession ", id)
    ind
  })
  if (any(duplicated(matches))) stop("RCC file(s) matching multiple accessions")
  matches
}


#' Extract the counts matrix from a `NanoTube::read_merge_rcc` result
#'
#' @param ns.data [NanoTube::read_merge_rcc()] result
#' @param colnames.func Function for generating new column names, given the "default" ones
get.counts.matrix <- function(ns.data, colnames.func=NULL) {
  counts <- ns.data$exprs
  rownames(counts) <- ns.data$dict[rownames(counts), "Name"]
  if (!is.null(colnames.func)) colnames(counts) <- colnames.func(colnames(counts))
  counts
}


plot.counts <- function(counts, row, coldata) {
  plot.data <- cbind(counts = counts[row, ], coldata)
  ## avoid "Don't know how to automatically pick scale ..." warning:
  plot.data$days <- as.numeric(plot.data$timepoint, units="days")
  gene <- ifelse(is.character(row), row, rownames(counts)[row])
  ggplot(plot.data, aes(days, counts, color = subject)) +
    geom_point(size = 2) +
    geom_line() +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5)) +
    labs(title = gene, subtitle = paste0("(", gene.categories[gene], ")"))
}

plot.all.counts <- function(counts, coldata, path) {
  stopifnot(all(colnames(counts) == rownames(coldata)))
  pdf(path)
  for (i in 1:nrow(counts)) {
    print(plot.counts(counts, i, coldata))
  }
  dev.off()
}

## PCA plot for control genes:
plotPCA.control <- function(counts, coldata, control.genes, title, legend.row=16,
                            color.batch=FALSE) {
  old.par <- par(mar=c(5, 4, 5, 2) + 0.1, xpd=TRUE)
  ## color by RCC run or batch?
  color.by <- if (color.batch) coldata$rcc_batch else coldata$rcc_run
  n.runs <- nlevels(color.by)
  colors <- scales::hue_pal()(n.runs)
  ## use only control genes?
  if (!is.null(control.genes)) {
    counts <- counts[control.genes, ]
    title <- paste0(title, ", control genes")
  }
  EDASeq::plotPCA(counts, labels=FALSE, col=colors[as.integer(color.by)], pch=19)
  title(title, line=3)
  ## legend parameters:
  single.row <- n.runs <= legend.row
  inset <- -0.8 / par("fin")[2]
  if (single.row) {
    what <- if (color.batch) "Batch:" else "Run:"
    legend("top", c(what, levels(color.by)), pch=c(NA, rep(19, n.runs)),
           col=c(NA, colors), bty="n", inset=c(0, inset), horiz=TRUE,
           ## use less horizontal space:
           x.intersp=0.5, text.width=NA)
  } else {
    n.rows <- ceiling(n.runs/legend.row)
    legend("top", levels(color.by), pch=19, col=colors, bty="n",
           ncol=ceiling(n.runs / n.rows), inset=c(0, inset - 0.01 * n.rows),
           ## use less horizontal space:
           x.intersp=0.5, text.width=NA)
  }
  par(old.par)
}

## RLE plots:
plotRLE <- function(counts, color.by=NULL, labels=NULL, range=c(-2, 2)) {
  plot <- ruv::ruv_rle(t(log2(counts + 1)), color.by, ylim = range)
  if (!is.null(labels))
    plot <- plot + geom_text(aes(label=labels), y=min(range) + 0.1, angle=90, hjust=0)
  plot
}

## MA plots:
plotMA.annotated <- function(res, dds, legend.pos="topleft", ...) {
  if (class(dds) == "DESeqDataSet") {
    ## check overlap
    stopifnot(all(sort(rownames(res)) == sort(rownames(rowData(dds)))))
    if (!(class %in% names(rowData(dds)))) stop("Missing 'class' column with gene class")
    ann <- rowData(dds)$class
    names(ann) <- rownames(rowData(data.ds))
  } else ann <- dds # expect vector classes named by genes

  DESeq2::plotMA(res, cex=1, ...)
  ## highlight housekeeping genes:
  ind <- ann[rownames(res)] == "Housekeeping"
  with(res[ind, ], points(baseMean, log2FoldChange, cex=2, lwd=3, col="orange"))
  ## highlight negative controls:
  ind <- ann[rownames(res)] == "Negative"
  with(res[ind, ], points(baseMean, log2FoldChange, cex=2, lwd=3, col="red"))
  ## highlight negative controls:
  ind <- ann[rownames(res)] == "Positive"
  with(res[ind, ], points(baseMean, log2FoldChange, cex=2, lwd=3, col="green"))

  legend(legend.pos, c("significant", "housekeeping", "neg. control", "pos. control"),
         col=c("blue", "orange", "red", "green"), pch=c(19, 1, 1, 1))
}


plot.top.genes <- function(plot.data, genes, qvalues, pvalues = NULL, n = 6, highlight=c(0, 0),
                           big.labels = TRUE) {
  ## q-values may be less granular than p-values, so sort by p-value if available:
  sel <- if (is.null(pvalues)) order(qvalues)[1:n] else order(pvalues)[1:n]
  hits <- genes[sel]
  qvalues <- qvalues[sel]
  ## some q-values may be "NA", e.g. due to "independent filtering"
  ## (https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#indfilt):
  labels <- ifelse(qvalues < 0.001, "q < 0.001", paste("q =", round(qvalues, 3)))
  labels <- paste0(hits, " (", labels, ")")
  names(labels) <- hits

  plot.part <- plot.data[plot.data$gene %in% hits, ]
  plot.part$gene <- factor(plot.part$gene, levels=hits) # preserve ordering
  plot.part$days <- as.numeric(plot.part$timepoint, units="days")
  baseline <- plot.part[plot.part$timepoint == 0, ]

  plot <- ggplot(plot.part, aes(days, counts, color=subject)) +
    geom_rect(aes(xmin=highlight[1], xmax=highlight[2], ymin=-Inf, ymax=Inf),
              color="transparent", fill="lightyellow") +
    geom_point(size = 2) +
    geom_line() +
    geom_hline(aes(yintercept=counts, color=subject), data=baseline, linetype="dashed") +
    facet_wrap(vars(gene), scales = "free_y", labeller=as_labeller(labels)) +
    theme_bw()

  if (big.labels) {
    plot + theme(strip.text=element_text(size=12, face="bold"))
  } else {
    plot + theme(strip.text=element_text(size=8))
  }
}


extract.terms1 <- function(str) {
  if (grepl("^cohort", str)) { # comparing cohorts
    unlist(strsplit(sub("^cohort(\\d)(\\.visit[ _]short(.*))?$", "\\1 \\3", str), " ", fixed=TRUE))
  } else { # comparing visits
    gsub(" ", "_", gsub("^visit short ", "", str))
  }
}

extract.terms <- function(res) {
  desc <- res@elementMetadata@listData$description[2]
  contrast <- sub("[^:]+: ", "", desc)
  if (grepl(" vs ", contrast)) { # comparison of two effects
    parts <- unlist(strsplit(contrast, " vs ", fixed=TRUE))
    terms <- lapply(parts, extract.terms1)
    stopifnot(length(terms[[1]]) == length(terms[[2]]))
    ## if timepoint (visit) is included, it should be the same on both sides:
    stopifnot((length(terms[[1]]) == 1) || (terms[[1]][2] == terms[[2]][2]))
    return(c(terms[[1]][1], terms[[2]])) # two cohorts, maybe with timepoint (visit), or two visits
  } else { # single effect
    terms <- extract.terms1(sub(" effect$", "", contrast))
    if (length(terms) == 1) return(c(terms, "1")) # cohort only
    return(c(terms[1], terms[1], terms[2])) # one cohort with timepoint (visit)
  }
}

## alternative version using ggplot2:
ggplot.top.genes <- function(res, counts, coldata, n.genes=12, use.time=FALSE, sort=TRUE,
                             show.legend=TRUE, focus.visit=NULL) {
  if (sort) res <- res[order(res$pvalue), ]
  genes <- rownames(res)[1:n.genes]
  labels <- paste0(genes, "\nLFC: ", format(res$log2FoldChange[1:n.genes], digits=2),
                   ", q: ", format(res$padj[1:n.genes], digits=3))
  names(labels) <- genes
  if (is.null(focus.visit)) {
    ## extract information on comparison from metadata:
    terms <- tryCatch(extract.terms(res), error=function(e) c("", ""))
    ind <- coldata$cohort == terms[1]
    if (!any(ind)) ind <- TRUE # no cohort selected -> use all data
    if (length(terms) > 2) { # visit-specific comparison
      focus.visit <- terms[3]
    } else if (isTRUE(ind)) { # comparison of time points, all cohorts
      focus.visit <- terms
    }
  } else ind <- TRUE
  if (length(focus.visit) == 1) focus.visit <- c("C1D1", focus.visit)
  plot.data <- lapply(genes, function(gene)
    data.frame(gene=gene, counts=counts[gene, ind],
               coldata[ind, c("timepoint", "subject", "cohort", "visit_short")]))
  plot.data <- do.call(rbind, plot.data)
  ## order genes by p-value in the plot, not by name:
  plot.data$gene <- factor(plot.data$gene, levels=genes)
  plot.part <- plot.data[plot.data$visit_short %in% focus.visit, ]
  xvar <- ifelse(use.time, "timepoint", "visit_short")
  plot <- ggplot(plot.data, aes(.data[[xvar]], counts)) +
    geom_line(aes(group=subject, color=subject), linetype=ifelse(nrow(plot.part) == 0, "solid", "dotted")) +
    geom_point(aes(color=subject), data=plot.part) +
    geom_line(aes(group=subject, color=subject), data=plot.part) +
    facet_wrap(vars(gene), scales="free_y", labeller=labeller(gene=labels)) +
    scale_y_log10() +
    theme(legend.position=ifelse(show.legend, "top", "none"), legend.key.spacing.y=unit(0, "cm")) +
    labs(title="Top differentially expressed genes", x=ifelse(use.time, "days", "visit"))
  if (!use.time) # rotate x axis labels for readability
    plot <- plot + theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))
  plot
}


expression.heatmap <- function(counts, gene.categories, sort=TRUE, log=TRUE, ...) {
  if (log) counts <- log10(counts + 1)
  if (sort) {
    totals <- rowSums(counts)
    counts <- counts[order(gene.categories, -totals), ]
    tab <- table(gene.categories)[sort(unique(gene.categories))]
    gaps <- cumsum(tab[1:(length(tab) - 1)])
  } else gaps <- NULL
  ## transpose matrix so dimensions fit better on a wide screen:
  pheatmap(t(counts), cluster_rows=FALSE, cluster_cols=FALSE, show_colnames=FALSE, gaps_col=gaps,
           annotation_col=data.frame(category=gene.categories[rownames(counts)]), ...)
}


#' Ratio-based scaling using standard (reference) counts
#'
#' @param counts Matrix of counts for samples of interest
#' @param std.counts Matrix of counts for standard samples
#' @param mapping Vector of indices mapping standard samples to the corresponding samples of interest
#' @param ref.counts Matrix of counts for reference RNA samples (optional)
ratio.scaling <- function(counts, std.counts, mapping, ref.counts=NULL) {
  stopifnot(all(std.counts > 0))
  stopifnot(ncol(counts) == length(mapping))
  stopifnot(all(unique(mapping) %in% 1:ncol(std.counts)))
  ## scale standard counts using reference RNA?
  if (!is.null(ref.counts)) {
    ref.genes <- rownames(ref.counts)
    ratios <- ref.counts / std.counts[ref.genes, ]
    medians <- apply(ratios, 1, median) # median ref. level per gene
    median.ratios <- apply(ratios, 2, function(col) col / medians)
    scale.factors <- apply(median.ratios, 2, median) # per run
    std.counts <- t(t(std.counts) * scale.factors)
  }

  medians <- apply(counts, 1, median) # median expression per gene
  ## scale each sample by its corresponding panel standard (reference sample):
  for (i in 1:ncol(counts)) {
    counts[, i] <- counts[, i] / std.counts[rownames(counts), mapping[i]]
  }
  ## scale each gene back to its original median (to keep magnitude for DGE analysis):
  new.medians <- apply(counts, 1, median)
  round(counts * medians / new.medians)
}


#' Generate model matrix for DESeq2
#'
#' Remove all-zero columns to avoid "model matrix is not full rank" error.
make.model.matrix <- function(formula, coldata) {
  model.mat <- model.matrix(formula, coldata)
  all.zero <- apply(model.mat, 2, function(col) all(col == 0))
  model.mat[, !all.zero]
}


#' Check whether DESeq2 model converged for all genes
check.convergence <- function(data.ds) {
  nonconv <- rownames(mcols(data.ds))[!mcols(data.ds)$betaConv]
  all.good <- length(nonconv) == 0
  if (!all.good) {
    message("DESeq2 failed to converge for: ", paste(nonconv, collapse=", "))
  } else message("No DESeq2 convergence failures")
  invisible(all.good)
}


#' Get gene sets defined by NanoString for IO 360 panel
#'
#' @param file Input Excel file with gene set information from NanoString
#' @param io.sheet Name of the sheet containing immuno-oncology-related annotations
#' @param func.sheet Name of the sheet containing functional annotations
#' @return List of gene sets in different categories
get.io360.gene.sets <- function(file="c:/Users/hendrik.weisser/OneDrive - Storm Therapeutics/Data/STC15_Nanostring/2024-09-19_New_LBL-10498-02_IO_360_Gene_List.xlsx",
                                io.sheet="Cancer-Immunity Cycle", func.sheet="Functional Annotations") {
  io.ann <- as.data.frame(read_excel(file, sheet=io.sheet, skip=1))
  func.ann <- as.data.frame(read_excel(file, sheet=func.sheet, skip=1))
  ## remove footers:
  io.ann <- io.ann[1:770, ]
  func.ann <- func.ann[1:770, ]
  ## cell types:
  part <- io.ann[!is.na(io.ann$"Cell Type"), c("Gene", "Cell Type")]
  ct.sets <- split(part$Gene, part$"Cell Type")
  ## "cancer-immunity cycle":
  ## TODO: this contains "Internal Reference Genes" as last set - remove (or keep as "negative control")?
  io.sets <- lapply(3:ncol(io.ann), function(col) io.ann$Gene[which(io.ann[[col]] == "+")])
  names(io.sets) <- names(io.ann)[-(1:2)]
  ## functional:
  func.sets <- lapply(3:ncol(func.ann), function(col) func.ann$Gene[which(func.ann[[col]] == "+")])
  names(func.sets) <- names(func.ann)[-(1:2)]
  ## remove entries from columns that don't contain gene sets (e.g. "Human Gene" in mouse data):
  io.sets <- io.sets[lengths(io.sets) > 0]
  func.sets <- func.sets[lengths(func.sets) > 0]
  list(cancerimmunity=io.sets, functional=func.sets, celltypes=ct.sets)
}


get.io360.gene.info <- function(file="c:/Users/hendrik.weisser/OneDrive - Storm Therapeutics/Data/STC15_Nanostring/2024-09-19_New_LBL-10498-02_IO_360_Gene_List.xlsx", sheet="Gene and Probe Details") {
  gene.info <- as.data.frame(read_excel(file, sheet, skip=1))
  gene.info <- gene.info[-nrow(gene.info), ] # remove copyright notice
  gene.info$Category <- "Endogenous"
  ind <- which(gene.info[[1]] == "Internal Reference Genes")
  gene.info <- gene.info[-ind, ]
  rownames(gene.info) <- NULL
  gene.info$Category[ind:nrow(gene.info)] <- "Housekeeping"
  stopifnot(nrow(gene.info) == 770)
  gene.info
}


gsea.nanostring <- function(dge.results, io360.sets=get.io360.gene.sets(), genes=NULL, pvalue.cutoff=0.1,
                            out.prefix="", titles=c(celltypes="Cell type", cancerimmunity="Cancer immunity",
                                                    functional="Functional annotations")) {
  gsea.res.all <- lapply(names(io360.sets), function(set.name) {
    gsea.res <- lapply(dge.results, function(res) {
      scores <- get.gsea.input(res, genes)
      gsea.custom(scores, io360.sets[[set.name]], pvalue.cutoff=pvalue.cutoff)
    })
    names(gsea.res) <- names(dge.results)
    if (out.prefix != "") {
      pdf(paste0(out.prefix, "_", set.name, ".pdf"))
      for (res.name in names(dge.results)) {
        print(dotplot.direction(gsea.res[[res.name]]) +
              labs(title=paste0(titles[set.name], ": ", res.name)) +
              theme(plot.title=element_text(hjust=1, face="bold")))
      }
      dev.off()
    }
    gsea.res
  })
  names(gsea.res.all) <- names(io360.sets)
  gsea.res.all
}


#' Visualize LFCs in each response group
boxplot.response.lfcs <- function(data.ds, genes, categories=levels(colData(data.ds)$response),
                                  title="", res.neg=NULL, res.pos=NULL, signif=NULL, show.n=TRUE) {
  ## DGE results per response group for each time point:
  if (is.null(res.neg)) {
    re <- paste0("^response", categories[1], "\\.visit_short")
    coefs.neg <- grep(re, resultsNames(data.ds), value=TRUE)
    res.neg <- lapply(coefs.neg, function(coef) as.data.frame(results(data.ds, list(coef))))
    names(res.neg) <- sub("^.*visit_short", "", coefs.neg)
  }
  if (is.null(res.pos)) {
    re <- paste0("^response", categories[2], "\\.visit_short")
    coefs.pos <- grep(re, resultsNames(data.ds), value=TRUE)
    res.pos <- lapply(coefs.pos, function(coef) as.data.frame(results(data.ds, list(coef))))
    names(res.pos) <- sub("^.*visit_short", "", coefs.pos)
  }
  ## stopifnot(all(names(res.pos) == names(res.neg)))

  lfcs.neg <- sapply(res.neg, function(res) res[genes, "log2FoldChange"])
  lfcs.pos <- sapply(res.pos, function(res) res[genes, "log2FoldChange"])

  lfcs.neg <- as.data.frame(tidyr::pivot_longer(data.frame(lfcs.neg, gene=genes), cols=!gene))
  lfcs.pos <- as.data.frame(tidyr::pivot_longer(data.frame(lfcs.pos, gene=genes), cols=!gene))

  plot.data <- rbind(lfcs.neg, lfcs.pos)
  plot.data$response <- rep(categories, c(nrow(lfcs.neg), nrow(lfcs.pos)))
  if (!is.null(signif)) { # mark some visits as significant
    ind <- plot.data$name %in% signif
    plot.data$name[ind] <- paste0(plot.data$name[ind], "*")
  }
  subtitle <- paste(length(genes), "genes:", paste(strwrap(paste(genes, collapse=", "), 60), collapse="\n"))

  plot <- ggplot(plot.data, aes(name, value)) +
    geom_hline(yintercept=0) +
    labs(x="Visit", y=bquote(log[2]~fold~change~(relative~to~C1D1)), title=title, caption=subtitle) +
    theme_bw() +
    theme(plot.caption=element_text(hjust=0.5)) # legend.position="top",
  if (show.n) {
    tab <- as.data.frame(table(colData(data.ds)[, c("visit_short", "response")]))
    tab <- tab[tab$visit_short != "C1D1", ]
    plot <- plot +
      geom_text(aes(visit_short, label=Freq, color=response, hjust=ifelse(response == "responder", -0.1, 1.1)),
                y=-Inf, data=tab, vjust=-0.2, show.legend=FALSE)
  }
  if (length(genes) > 1) plot + geom_boxplot(aes(fill=response), position=position_dodge2(preserve="single"))
  else plot + geom_point(aes(color=response), size=3) + geom_line(aes(color=response, group=response))
}


get.enriched.genes <- function(description, gsea.results) {
  lapply(gsea.results, function(res) {
    if (nrow(res) == 0) return(NULL)
    genes <- res$core_enrichment[res$Description == description]
    sort(unlist(strsplit(genes, "/", fixed=TRUE)))
  })
}

common.enriched <- function(enriched.genes) {
  all.genes <- unique(unlist(enriched.genes))
  filtered <- enriched.genes[!sapply(enriched.genes, is.null)]
  common <- Reduce(intersect, filtered)
  message(length(common), "/", length(all.genes), " genes in common")
  common
}

get.pathway.genes <- function(collection, term=NULL, description=NULL, species="Homo sapiens") {
  ## need either 'term' or 'description' of gene set:
  stopifnot(!is.null(term) || !is.null(description))
  stopifnot(is.null(term) || is.null(description))
  if (toupper(collection) == "REACTOME") {
    if (is.null(term)) {
      map <- as.list(reactome.db::reactomePATHNAME2ID)
      ## Reactome(DB) adds species prefix (e.g. "Homo sapiens: "), so use 'grep':
      hits <- grep(description, names(map), fixed=TRUE, value=TRUE)
      if (length(hits) == 1) {
        term <- map[[hits]]
      } else if (length(hits) > 1) {
        species.desc <- paste0(species, ": ", description)
        if (species.desc %in% hits) {
          term <- map[[species.desc]]
        } else stop("Multiple matches in Reactome for: ", description)
      } else stop("No match in Reactome for: ", description)
    }
    gene.ids <- unlist(AnnotationDbi::mapIds(reactome.db::reactome.db, keys=term, column="ENTREZID",
                                             keytype="PATHID", multiVals="list"))
    gene.map <- clusterProfiler::bitr(gene.ids, "ENTREZID", "SYMBOL", org.Hs.eg.db::org.Hs.eg.db, FALSE)
    unique(na.omit(gene.map$SYMBOL))
  } else if (startsWith(collection, "GO")) {
    if (is.null(term)) term <- AnnotationDbi::mapIds(GO.db::GO.db, keys=description, column="GOID", keytype="TERM")
    gene.map <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys=term, columns=c("SYMBOL", "GENETYPE"),
                                      keytype="GOALL")
    unique(na.omit(gene.map[gene.map$GENETYPE == "protein-coding", "SYMBOL"]))
  }
}


capitalize.first <- function(str) {
  paste0(toupper(substr(str, 1, 1)), substr(str, 2, nchar(str)))
}


## for ambiguous gene probes, select the gene with more GO annotations:
fix.nanostring.genes <- function(genes, replace=c("CCL3/L1"="CCL3", "FCGR3A/B"="FCGR3A", "MAGEA3/A6"="MAGEA6",
                                                  "TPSAB1/B2"="TPSAB1", "XCL1/2"="XCL1"), fix.names=FALSE) {
  if (fix.names) {
    names(genes)[match(names(replace), names(genes))] <- replace
  } else {
    genes[match(names(replace), genes)] <- replace
  }
  genes
}


plot.qc.medians <- function(counts, pdf.path, gene.categories, subtitles=NULL, qc.failed=c()) {
  medians <- rowMedians(counts[, !(colnames(counts) %in% qc.failed)])
  pdf(pdf.path)
  on.exit(dev.off())
  for (sample in colnames(counts)) {
    title <- paste0(sample, " (QC ", ifelse(sample %in% qc.failed, "FAILED", "passed"), ")")
    cat.colors <- c(Housekeeping="yellow", Negative="red", Positive="green")
    plot(medians + 1, counts[, sample] + 1, log="xy", xlim=c(1, max(counts)), ylim=c(1, max(counts)),
         main=title, xlab="Median counts (QC-passed samples)", ylab="Sample counts")
    grid()
    abline(coef=c(0, 1))
    for (category in names(cat.colors)) {
      points(medians[gene.categories == category] + 1, counts[gene.categories == category, sample] + 1, pch=19,
             cex=0.7, col=cat.colors[category])
    }
    legend("topleft", names(cat.colors), pch=19, col=cat.colors, title="Controls:", bg="white")
    if (!is.null(subtitles)) mtext(subtitles[sample], line=0.5)
  }
}


## remove batch effects estimated using RUVg from a matrix of counts
remove.batch.effects <- function(counts, coldata, is.log=FALSE, min.zero=TRUE) {
  w.vars <- grep("^W_\\d+$", colnames(coldata))
  ## 'limma::removeBatchEffect' expects log-expression values:
  if (!is.log) counts <- log2(counts + 1)
  norm.counts <- limma::removeBatchEffect(counts, covariates=coldata[, w.vars])
  if (!is.log) norm.counts <- 2^norm.counts - 1 # undo log-transform and pseudocount
  if (min.zero) norm.counts[norm.counts < 0] <- 0
  norm.counts
}
