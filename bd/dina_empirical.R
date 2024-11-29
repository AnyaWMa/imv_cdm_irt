source("00funs.R")

#########################irw
dataset <- redivis::organization("datapages")$dataset("Item Response Warehouse")
# List all tables in the dataset
dataset_tables <- dataset$list_tables()
# Rename the tables using their names for easier access
names(dataset_tables) <- sapply(dataset_tables, function(x) x$name)

f<-function(tab) {
    df <- tab$to_data_frame()  # Convert the table to a data frame
    df$resp <- as.numeric(df$resp)  # Ensure the response variable is numeric
    resp<-irw::long2resp(df)
    resp$id<-NULL
    makeqm<-function(df,resp=NULL) { #resp required for correct ordering
        ii<-grep("Qmatrix__",names(df))
        L<-split(df,df$item)
        qm<-lapply(L,function(x,ii) colMeans(x[,ii],na.rm=TRUE),ii)
        qm<-do.call("rbind",qm)
        if (!is.null(resp)) {
            cn<-colnames(resp)
            qm<-qm[cn,]
        }
        qm
    }
    qm<-makeqm(df,resp=resp)
    list(resp=resp,qm=qm)
}

om<-list()
for (tables in c("frac20","cdm_ecpe","cdm_hr","mcmi_mokken")) {
    L<-f(dataset_tables[[tables]])
    om[[tables]]<-oos.compare(L$resp,L$qm,5,modeltype="2PL")
}
unlist(om)
