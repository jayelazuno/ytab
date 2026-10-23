#!/usr/bin/env Rscript
# Analysis scripts write no progress events, so the job manager closes their record on exit: every queued item takes
# the exit status, counters are complete, and the elapsed time comes from the recorded start (not 0/N "queued").
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1L)
`%||%`<-function(x,y)if(is.null(x)||!length(x))y else x
repo_root<-normalizePath(args[[1L]],mustWork=TRUE)
source(file.path(repo_root,"app/shiny/R/job_manager.R"))
p<-list(status="running",total_items=2L,processed_items=0L,successful_items=0L,skipped_items=0L,failed_items=0L,remaining_items=2L,
        job_started_at="2026-10-21 04:13:52 UTC",items=list(list(item="betageom",status="queued"),list(item="ttr",status="queued")))
now<-as.POSIXct("2026-10-21 04:15:31",tz="UTC")
ok<-job_progress_close(p,"success",now);bad<-job_progress_close(p,"failed",now)
stopifnot(ok$status=="success",ok$processed_items==2,ok$successful_items==2,ok$remaining_items==0,ok$job_elapsed_seconds==99,
          all(vapply(ok$items,`[[`,"","status")=="success"),ok$message=="Finished.",
          bad$failed_items==2,bad$successful_items==0,all(vapply(bad$items,`[[`,"","status")=="failed"))
legacy<-c(p,list(status="success",message="Process ended before a final progress event.",updated_at="2026-10-21 04:15:31 UTC"));legacy$status<-"success"
shown<-job_progress_display(legacy);stopifnot(shown$processed_items==2,shown$job_elapsed_seconds==99,identical(job_progress_display(ok),ok))
cat("PASS: progress records of scripts without progress events close with complete counters and real elapsed time\n")
