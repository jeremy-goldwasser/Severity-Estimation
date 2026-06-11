library(batchtools)
loadRegistry("real-data-reg")
res <- reduceResultsList()
res$meta <- getJobTable() |> unwrap()

write_rds(res, here::here("Data", "batch-results", "real-data-res.rds"))
