options(warn =-1)
options(tidyverse.quiet = TRUE)
library("tidyverse")

Args = commandArgs()
infile = Args[6]
sample=Args[7]
recipientName=Args[8]
donorName=Args[9]
tab = read.table(infile,header=T,sep="\t",stringsAsFactors = F)

#======pre-process======
tab=tab[tab$All.dp>100,]
tab=tab[tab$Chr!="chrX",]
tab=tab[tab$Chr!="chrY",]
print(paste0("after removing chrxy, and dp<500:",dim(tab)[1]))
tab <- tab %>% filter((Donor.rf < quantile(Donor.rf)[4] + 1.5 * IQR(Donor.rf)) & (Donor.rf > quantile(Donor.rf)[2]-1.5 * IQR(Donor.rf)),na.rm = TRUE)
print(paste0("after removing outlier:",dim(tab)[1]))                   

rf.min = 0.00125


#======remove outlier by MAD========
smoothcount<-function(hcount)
{
  n = length(hcount)
  hcount2 = hcount[1]
  barplot(height = hcount,width = 2,col = "black")
  for(i in 2:(n-1))
  {
    m = round(mean(hcount[(i-1):(i+1)]),3)
    hcount2 = c(hcount2,m) 
  }
  hcount2 = c(hcount2,hcount[n])
  barplot(height = hcount2,width = 2,col = "red")
  return(hcount2)
}

hist_range_for_median_calculation <- function(hcount,min,max)
{
  sums = cumsum(hcount)
  start = which(sums>=min)[1]
  temp = which(sums<=max)
  end = temp[length(temp)] + 1
  return(c(start,end))
}

newmedian <- function(nn)
{
	h = hist(nn,breaks=50,plot=FALSE)
	rf = h$breaks
	hcount = smoothcount(h$count)
	n = ifelse(length(rf)>length(hcount),length(hcount),length(rf))
	#min = round(length(nn)/4*1,0)
	#max = round(length(nn)/4*3,0)
	min = length(nn[nn<=0])
	max = round(length(nn)/4*3,0)
	range = hist_range_for_median_calculation(hcount,min,max)
	#print(c(n,start,end))
	rf.range = rf[range[1]:range[2]]
	hcount.range = hcount[range[1]:range[2]]
	mix = rf.range * hcount.range
	c=data.frame(rf=rf.range,hcount=hcount.range,mix=mix)
	
	median_classic = median(nn)
	median_maxrf = rf.range[mix==max(mix)]
	median_maxhcount = rf.range[hcount.range==max(hcount.range[rf.range>rf.min])]
	#median = ifelse(median_classic<=median_new, median_new, median_classic)
	median = max(c(median_classic,median_maxrf[1],median_maxhcount[length(median_maxhcount)]))
	
	
	h2=hist(nn,breaks=100,plot=FALSE)
	rf = h2$breaks
	hcount = smoothcount(h2$count)
	n2 = ifelse(length(rf)>length(hcount),length(hcount),length(rf))
	range = hist_range_for_median_calculation(hcount,min,max)
	rf.range = rf[range[1]:range[2]]
	hcount.range = hcount[range[1]:range[2]]
	mix = rf.range * hcount.range
	c2=data.frame(rf=rf.range,hcount=hcount.range,mix=mix)
	
	maxvalue_of_mix = max(c2[c2$rf>median*0.7&c2$rf<median*1.3,]$mix)
	median_maxhrf = c2[c2$mix==maxvalue_of_mix,]$rf
	median_maxhcount = rf.range[hcount.range==max(hcount.range[rf.range>rf.min])
	median = max(c(median,median_maxhrf[1],median_maxhcount[length(median_maxhcount)]))
	
	hist(nn,breaks=100)
	abline(v=median_classic,col="black")
	abline(v=median_maxhrf,col="red")
	abline(v=median_maxhcount,col="blue")


	return(median)
}

mad <-function(tab,n,nmedian)
{
	#tab: *mp.tab.ex
	#n: multiplicator or MAD

	deviations = abs(tab$Donor.rf - nmedian)
	mad = median(deviations)
	#tag = rep("drop",length(tab$Donor.rf))
	#tag[(tab$Donor.rf-nmedian)<n*mad] = "keep"
	#print(paste0("rf: ",tab$Donor.rf-median, " median: ",median," deviations: ",deviations, " mad: ", mad, " n*mad: ", n*mad," ",tag))
	newd = tab[deviations<n*mad,]
	return(newd)
}

nmedian = newmedian(tab$Donor.rf-tab$Noise.rf)

tab = mad(tab,2,nmedian)


#=====removing >max & <min============
tab=tab[order(tab$Donor.rf,decreasing=F),]
n = length(tab$Donor.rf)


#======initial ddcfDNA=====
nmedian = newmedian(tab$Donor.rf)

rf.max = ifelse(nmedian * 3 < 1, nmedian * 3, 1)

#======population MAF======
popF = data.frame(Homo=tab$D.homo.AF,Het=tab$D.hete.AF)


p.ddcfDNA.by_twogenome_multiplexPCR <- function(df)
{
	donordp = df$D.hete.AF * df$Donor.dp * 2 + df$D.homo.AF * df$Donor.dp
	totaldp = df$Recipient.dp+ df$Donor.dp
	donor.rfs = donordp / totaldp
	donor.rfs = donor.rfs[donor.rfs>rf.min]
	rf.max = mean(donor.rfs)*5
	donor.rfs = round(donor.rfs[donor.rfs<rf.max],3)
	
	donor.rfdf = data.frame(table(donor.rfs))
	donor.rfdf = donor.rfdf[order(donor.rfdf$Freq,decreasing=TRUE),]
	donor.percent = as.numeric(as.character(donor.rfdf[1,1]))
	return(donor.percent)
}

#======ddcfDNA.by_maxlike======
p.ddcfDNA.by_maxlike <- function(df,popdf){
	samplename = as.character(df[1,1])
	
	p.min = rf.min
	p.max = rf.max
	e.max = 0.005
	e.max = ifelse(e.max>nmedian,nmedian,e.max)

	k.range = c(p.min, p.max, nmedian)
	e.mu = c(-0.001, 0.001, -0.000247562)
	e.sigma = c(0, 0.0006, 0.0002)
	e.range = c(0.00125,e.max,0.0025)

	k.range[1] = k.range[3]
	k.range[2] = nmedian*3 + 0.001
	k.range[1] = ifelse(k.range[1] < p.min, p.min, k.range[1])
	k.range[2] = ifelse(k.range[2] > p.max, p.max, k.range[2])
	write(file=paste0(samplename,".likehood.log"),x=c("ddcfDNA.percent","p.e","p.homo","p.het","likehood"),sep="\t",append=FALSE)

	#probability model
	prob_model = function(para)
	{
		#=============binomial distribution====================
		d.Aa = dbinom(x=df$Donor.dp, size=(df$Recipient.dp+df$Donor.dp), prob=para[1]/2)
		d.AA = dbinom(x=df$Donor.dp, size=(df$Recipient.dp+df$Donor.dp), prob=para[1])
		d.ee = dbinom(x=df$Donor.dp, size=(df$Recipient.dp+df$Donor.dp), prob=para[4])

		#d.ee = dnorm(x=df$Donor.rf-df$Noise.rf, mean=para[2], sd=para[3]) 
		#d.ee = dnorm(x=df$Donor.rf, mean=-0.000247562,sd=0.00057102)
		
		p.e =  d.ee * df$Noise.rf
		p.het = popdf$Het*d.Aa
		p.homo = popdf$Homo*d.AA
		p.all = p.e + (p.het + p.homo) * (1-df$Noise.rf)

		max.like=sum(log(p.all[p.all > 0]))

		return(-max.like)
	}

	#esitimation
	model.est = nlminb(
		start=c(k.range[3],e.mu[3],e.sigma[3],e.range[3]),#initial value
		objective=prob_model,#probability model
		lower = c(k.range[1],e.mu[1],e.sigma[1],e.range[1]),#lower
		upper = c(k.range[2],e.sigma[2],e.sigma[2],e.range[2]),#upper
		)
	

	#related or unrelated transplant
	initialvalue = k.range[3]
	prediction = round(model.est$par[1],4)

	Rindex = prediction / initialvalue
	if(Rindex > 1.5)
	{
		print("prediction of relationship: Unrelated")
	}else{
		print("prediction of relationship: Related")
		prediction = prediction * 2
	}

	# return
	data.frame(sample = samplename,#
				recipient = recipientName,#
				allelscount = nrow(df), #records of SNP for percent estimation
				initial.value = k.range[3],
				ddcfDNA.percent = prediction)
}

#ddcfDNA prediction
result.df=p.ddcfDNA.by_maxlike(df=tab,popdf=popF)
write.table(x=result.df,file=paste0(sample,".ddcfDNA_percent.xls"),quote=F,sep="\t",row.names = F)