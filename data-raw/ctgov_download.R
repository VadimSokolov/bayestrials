# Port of code/ctgov_stats.py download step (sections 1 and 2) using base R and jsonlite.
# Requires jsonlite (not a package dependency). Run from the book root.
# Rscript bayestrials/data-raw/ctgov_download.R [--refresh]; then run ctgov.R.
# Note: the Python script also builds the by-area/phase/sponsor tables and figures; the
# tables are one-line table() calls on ctgov_bayesian_records and figures are drawn in the book.
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("install jsonlite to run the download")
base <- "https://clinicaltrials.gov/api/v2/studies"
D <- "data/ctgov"
dir.create(D, showWarnings = FALSE, recursive = TRUE)
refresh <- "--refresh" %in% commandArgs(TRUE)
today <- format(Sys.Date())
terms <- c(
  "Bayesian" = '"Bayesian"',
  "Adaptive design/trial" = '"adaptive design" OR "adaptive trial"',
  "Response-adaptive/adaptive randomization" = '"response-adaptive" OR "adaptive randomization"',
  "Group sequential/interim analysis" = '"group sequential" OR "interim analysis"',
  "Platform trial/master protocol" = '"platform trial" OR "master protocol"',
  "Predictive probability" = '"predictive probability"',
  "Thompson sampling" = '"Thompson sampling"')
interv <- "AREA[StudyType]INTERVENTIONAL"
ph23 <- "AREA[Phase](PHASE2 OR PHASE3)"
years <- 2000:2025

get_json <- function(params) {
  Sys.sleep(0.2)
  url <- paste0(base, "?", paste(names(params), vapply(params, function(v) utils::URLencode(as.character(v), reserved = TRUE), ""), sep = "=", collapse = "&"))
  for (k in 1:5) {
    out <- tryCatch(jsonlite::fromJSON(url, simplifyVector = FALSE), error = function(e) NULL)
    if (!is.null(out)) return(out)
    Sys.sleep(2 * k)
  }
  stop("request failed: ", url)
}
count <- function(term, design, year) {
  adv <- if (design == "phase23") paste(interv, "AND", ph23) else interv
  adv <- sprintf("%s AND AREA[StartDate]RANGE[%d-01-01,%d-12-31]", adv, year, year)
  p <- list(filter.advanced = adv, countTotal = "true", pageSize = 1, fields = "NCTId")
  if (!is.null(term)) p$query.term <- term
  get_json(p)$totalCount
}

cf <- file.path(D, "raw_counts.json")
if (file.exists(cf) && !refresh) {
  cache <- jsonlite::fromJSON(cf, simplifyVector = FALSE)
} else {
  cache <- list(download_date = today, counts = list())
  for (design in c("interventional", "phase23")) {
    for (nm in c("ALL", names(terms))) {
      term <- if (nm == "ALL") NULL else terms[[nm]]
      for (y in years) cache$counts[[paste(design, nm, y, sep = "|")]] <- count(term, design, y)
      message(design, " ", nm, " done")
    }
  }
  writeLines(jsonlite::toJSON(cache, auto_unbox = TRUE, pretty = TRUE), cf)
}
key <- strsplit(names(cache$counts), "|", fixed = TRUE)
df <- data.frame(design = vapply(key, `[`, "", 1), term = vapply(key, `[`, "", 2),
                 year = as.integer(vapply(key, `[`, "", 3)), n = unlist(cache$counts), stringsAsFactors = FALSE)
den <- df[df$term == "ALL", c("design", "year", "n")]
names(den)[3] <- "denominator"
cnt <- merge(df[df$term != "ALL", ], den, by = c("design", "year"))
cnt$share_pct <- 100 * cnt$n / cnt$denominator
cnt$download_date <- cache$download_date
cnt <- cnt[order(match(cnt$design, c("interventional", "phase23")), match(cnt$term, names(terms)), cnt$year),
           c("design", "term", "year", "n", "denominator", "share_pct", "download_date")]
write.csv(cnt, file.path(D, "counts_by_year.csv"), row.names = FALSE)

rf <- file.path(D, "raw_bayesian_records.json")
if (file.exists(rf) && !refresh) {
  rec <- jsonlite::fromJSON(rf, simplifyVector = FALSE)
} else {
  fields <- paste0("NCTId,BriefTitle,Phase,LeadSponsorClass,Condition,StartDate,OverallStatus,",
                   "EnrollmentCount,BriefSummary,DetailedDescription")
  studies <- list(); tok <- NULL
  repeat {
    p <- list(query.term = terms[["Bayesian"]], filter.advanced = interv, pageSize = 500,
              fields = fields, countTotal = "true")
    if (!is.null(tok)) p$pageToken <- tok
    j <- get_json(p)
    studies <- c(studies, j$studies)
    tok <- j$nextPageToken
    if (is.null(tok)) break
  }
  rec <- list(download_date = today, query = terms[["Bayesian"]], studies = studies)
  writeLines(jsonlite::toJSON(rec, auto_unbox = TRUE), rf)
}
areas <- list(
  Oncology = "cancer|carcinoma|tumou?r|neoplasm|leukemi|lymphoma|myeloma|melanoma|sarcoma|glioma|glioblastoma|metasta|oncolog|malignan",
  Cardiovascular = "cardi|heart|coronary|myocard|hypertens|stroke|atrial|vascular|thromb|embol|infarct|aort",
  "Infectious disease" = "infect|covid|sars|hiv|virus|viral|bacter|sepsis|pneumonia|tubercul|malaria|hepatitis|influenza|pathogen|vaccin",
  "Neurology/psychiatry" = "alzheimer|parkinson|epilep|neuro|dementia|depress|schizo|psych|migraine|multiple sclerosis|brain|anxiety|cognitive|pain")
area_of <- function(conds) {
  t <- tolower(paste(conds, collapse = " "))
  for (a in names(areas)) if (grepl(areas[[a]], t, perl = TRUE)) return(a)
  "Other"
}
chr <- function(x, default = "") if (is.null(x)) default else x
out <- lapply(rec$studies, function(s) {
  pm <- s$protocolSection
  conds <- unlist(pm$conditionsModule$conditions)
  sd <- chr(pm$statusModule$startDateStruct$date)
  ph <- paste(unlist(pm$designModule$phases), collapse = "/")
  data.frame(nct = pm$identificationModule$nctId, title = chr(pm$identificationModule$briefTitle),
             phase = if (ph == "") "NA" else ph,
             sponsor_class = chr(pm$sponsorCollaboratorsModule$leadSponsor$`class`, "UNKNOWN"),
             conditions = paste(conds, collapse = "; "), area = area_of(conds), start_date = sd,
             start_year = if (nzchar(sd)) as.numeric(substr(sd, 1, 4)) else NA_real_,
             status = chr(pm$statusModule$overallStatus),
             enrollment = chr(pm$designModule$enrollmentInfo$count, NA_real_), stringsAsFactors = FALSE)
})
write.csv(do.call(rbind, out), file.path(D, "bayesian_records.csv"), row.names = FALSE)
