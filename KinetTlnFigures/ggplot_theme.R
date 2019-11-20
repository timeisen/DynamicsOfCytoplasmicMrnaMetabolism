#Tim's theme
#Created 20171129

theme_tim <- function(base_size = 6, base_family = "Helvetica")
{
theme_classic(base_size = base_size, base_family = base_family) +
theme(
	panel.grid=element_blank(),
	axis.text.y = element_text (size = 6,color="black"),
	axis.text.x = element_text (size = 6,color="black"),
	axis.ticks.length = unit(0.1 , "cm"),
	#axis.line=element_line(size=.5,color="black"),
	#axis.text.x=element_blank(),
  #axis.ticks=element_blank(),
  #axis.title.x=element_blank(),
  legend.position="none",
  strip.background = element_blank(),
  panel.spacing=unit(0.3,"cm"),
  # axis.line=element_blank(),
  # axis.text.x=element_blank(),
  # axis.ticks.x=element_blank(),
  # axis.text.y=element_blank(),
  axis.title.y=element_blank(),
  axis.title.x=element_blank(),
  strip.text.x = element_text(size = 6, color = "black"))
}

theme_tim_label <- function(base_size = 6, base_family = "Helvetica")
{
theme_classic(base_size = base_size, base_family = base_family) +
theme(
  panel.grid=element_blank(),
  axis.text.y = element_text (size = 6,color="black"),
  axis.text.x = element_text (size = 6,color="black"),
  axis.ticks.length = unit(0.1 , "cm"),
  #axis.line=element_line(size=.5,color="black"),
  #axis.text.x=element_blank(),
  #axis.ticks=element_blank(),
  #axis.title.x=element_blank(),
  # legend.position="none",
  strip.background = element_blank(),
  panel.spacing=unit(0.3,"cm"),
  # axis.line=element_blank(),
  # axis.text.x=element_blank(),
  # axis.ticks.x=element_blank(),
  # axis.text.y=element_blank(),
  # axis.title.y=element_blank(),
  # axis.title.x=element_blank(),
  strip.text.x = element_text(size = 6, color = "black"))
}

#Larger font, for presentations
theme_tim_presentation <- function(base_size = 12, base_family = "Helvetica")
{
theme_classic(base_size = base_size, base_family = base_family) +
theme(
  panel.grid=element_blank(),
  axis.text.y = element_text (size = 12,color="black"),
  axis.text.x = element_text (size = 12,color="black"),
  axis.ticks.length = unit(0.3 , "cm"),
  #axis.line=element_line(size=.5,color="black"),
  #axis.text.x=element_blank(),
  #axis.ticks=element_blank(),
  #axis.title.x=element_blank(),
  # legend.position="none",
  strip.background = element_blank(),
  panel.spacing=unit(0.3,"cm"),
  # axis.line=element_blank(),
  # axis.text.x=element_blank(),
  # axis.ticks.x=element_blank(),
  # axis.text.y=element_blank(),
  # axis.title.y=element_blank(),
  # axis.title.x=element_blank(),
  strip.text.x = element_text(size = 12, color = "black"))
}
# log_ticks<-function(start,end){
#   tic_num<-log10(end)-log10(start)+1
#   labels<-rep(" ",10*tic_num)
#   labels[seq(10,(tic_num*10),10)]<-10^(log10(start):log10(end))*10
#   breaks<-c(matrix(1:10)%*%t(matrix(10^(log10(start):log10(end)))))
#   return(list(labels,breaks))
# }


log_ticks<-function(start,end){
  tic_num<-log10(end)-log10(start)+1
  labels<-rep(" ",10*tic_num)
  labels[seq(10,(tic_num*10),10)]<-10^(log10(start):log10(end))*10
  breaks<-c(matrix(1:10)%*%t(matrix(10^(log10(start):log10(end)))))
  return(list(labels,breaks))
}

# log2_ticks<-function(start,end){ #What am I trying to do here
#   tic_num<-log2(end)-log2(start)+1
#   labels<-rep(" ",2*tic_num)
#   labels[seq(2,(tic_num*2),2)]<-2^(log2(start):log2(end))*2
#   breaks<-c(matrix(1:2)%*%t(matrix(2^(log2(start):log2(end)))))
#   return(list(labels,breaks))
# }

log_ticks_exp<-function(start,end){
  tic_num<-log10(end)-log10(start)+1
  labels<-rep(" ",10*tic_num)
  labels[seq(10,(tic_num*10),10)]<-(log10(start)+1):(log10(end)+1)
  breaks<-c(matrix(1:10)%*%t(matrix(10^(log10(start):log10(end)))))
  return(list(labels,breaks))
}

fancy_scientific <- function(l) {
     # turn in to character string in scientific notation
     l <- format(l, scientific = TRUE)
     #Don't display anything for minor ticks
     l <- gsub("2e.*","NA",l)
     l <- gsub("3e.*","NA",l)
     l <- gsub("4e.*","NA",l)
     l <- gsub("5e.*","NA",l)
     l <- gsub("6e.*","NA",l)
     l <- gsub("7e.*","NA",l)
     l <- gsub("8e.*","NA",l)
     l <- gsub("9e.*","NA",l)

     #display 0 as 0
     l <- gsub("0e\\+00","0",l)
     l <- gsub("1e","e",l)

     # quote the part before the exponent to keep all the digits
     l <- gsub("^(.*)e", "'\\1'e", l)
     # turn the 'e+' into plotmath format
     # remove + after exponent, if exists. E.g.: (3x10^+2 -> 3x10^2) 
     l <- gsub("e\\+","e",l) 
     # convert 1x10^ or 1.000x10^ -> 10^ 

     l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)

     l <- gsub("e", "%*%10^", l)
     #remove the X before the 10^
     l <- gsub("%*%", " ", l)
     # return this as an expression
     parse(text=l)
}

#Saving shortcut
ggsaveTim <- function(p){
  filename = paste0(deparse(substitute(p)),".pdf")
  print(filename)
  ggsave(plot = p, file = filename, width = 3, height = 3, useDingbats = FALSE)
}

cat(paste("log_ticks usage:","labels = log_ticks(start,end)[[1]]", "breaks = log_ticks(start,end)[[2]] \n",sep="\n"))