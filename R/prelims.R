#saving the data for use in "examples"

# st_code = "SCH"
# st = which(colnames(Rain)==st_code)
# sea = 2
# seas = seasons[[sea]]
#
# #aggretate the data
# precipdata =  seas[,c("date", "SCH")]
# save(precipdata, file = "data/precipdata.RData", compress = "xz")
#
#
# library(devtools)
# use_vignette(name = "eppdIDF_models",title =  "Modeling IDF Curves using the Extended Genaralized Pareto Distribution" )


#saving the data for use in "examples"

# st_code = "SCH"
# st = which(colnames(Rain)==st_code)
# sea = 2
# seas = seasons[[sea]]
#
# #aggretate the data
# precipdata =  seas[,c("date", "SCH")]
# save(precipdata, file = "data/precipdata.RData", compress = "xz")
#
#
# library(devtools)
# use_vignette(name = "eppdIDF_models",title =  "Modeling IDF Curves using the Extended Genaralized Pareto Distribution" )

##IDAF

# st_code ="BUS"
# sea = 1
#
# seas =  map_dfc(station_data_cpc_by_area, function(x) pluck(x,c(st_code))) %>%
#   mutate(date = station_data_cpc_by_area$`area_ 1`$date) %>% relocate(date) %>%
#   get_seasonal_data(sea)
#areal_data = seas
#save(areal_data, file = "data/precip_areal.RData", compress = "xz")
