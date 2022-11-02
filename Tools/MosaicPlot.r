library(ggplot2)



data=read.csv('SuperpopulationChrAll.PC20SVMResults',sep='\t')
metadata=read.csv('metadata.txt',sep='\t')
data=merge(metadata,data,by="run_accession")



ra_ids=c() 
pop=c()
prob=c()
self_reported=c()
for(i in 1:nrow(data)) {       # for-loop over rows
  ra_ids=append(ra_ids, replicate(5, data$run_accession[i]))
  pop=append(pop, c("AFR","AMR","EAS","EUR","SAS"))
  prob=append(prob,as.numeric(subset(data[rownames(data) %in% c(i),], select=AFR:SAS)))
  self_reported=append(self_reported, replicate(5, data$Eth1[i]))
}

df = data.frame(run_accession = ra_ids,
                Super.Population = pop, Self.Reported=self_reported,
                p = prob)

#Doesn't make diff(these next two lines) 
df$run_accession <- factor(df$run_accession, levels=unique(df$run_accession))
df$Super.Population <- factor(df$Super.Population)

df$Self.Reported=gsub('EAS', 'East Asian', df$Self.Reported)
df$Self.Reported=gsub('SAS', 'South Asian', df$Self.Reported)
df$Self.Reported=gsub('EUR', 'European', df$Self.Reported)
df$Self.Reported=gsub('AMR', 'American Admix', df$Self.Reported)
df$Self.Reported=gsub('AFR', 'African', df$Self.Reported)

ggplot(df, aes(x=run_accession, y=p, fill=Super.Population)) +
  geom_bar(stat="identity", position="stack",width = 1)   + xlab("RNA-Seq Samples") + ylab("RIA Ancestry Probability") + facet_grid(.~  Self.Reported, scales="free_x", space = "free_x") +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values=c("red",
                             "green",
                             "blue",
                             "purple","orange"))

ggsave("ADMIXTURE.png")



