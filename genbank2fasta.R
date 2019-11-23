seq<-readLines("AB115403.gb")
#1
a_n<-unlist(strsplit(seq[1],split = "\\s+"))[2]
a_n<-paste0(">",a_n)
#2
st<-grep("ORIGIN",seq)
ed<-grep("^//",seq)
seq1<-seq[(st+1):(ed-1)]
seq2<-paste(seq1,collapse = "")
seq3<-gsub("\\d","",seq2)
gb_seq<-gsub(" ","",seq3)
#3
fastafile<-c(a_n,gb_seq)
res<-paste0(a_n,"\n",gb_seq)
writeLines(res)
write.table(res,
            file = "test.fas",
            row.names = F,
            col.names = F,
            quote = F)
#test
readLines("test.fas")
library(seqinr)
read.fasta("test.fas")



#function1
read_genbank<-function(file){
  seq<-readLines(file)
  a_n<-unlist(strsplit(seq[1],split = "\\s+"))[2]
  st<-grep("ORIGIN",seq)
  ed<-grep("^//",seq)
  gb_seq<-gsub(" ","",gsub("\\d","",paste(seq[(st+1):(ed-1)],collapse = "")))
  count<-nchar(gb_seq)
  writeLines(paste0("accession number:","\n",a_n,"\n","\n","sequence:","\n",gb_seq,"\n","\n","sequence length:","\n",count,"\n","\n","Base compositon:","\n"))
  s<-tolower(unlist(strsplit(gb_seq,"")))
  count_A<-0
  count_T<-0
  count_G<-0
  count_C<-0
  for (j in 1:length(s)){
      if(s[j]=='a'){count_A=count_A+1}
      if(s[j]=='t'){count_T=count_T+1}
      if(s[j]=='g'){count_G=count_G+1}
      if(s[j]=='c'){count_C=count_C+1}
  }
  m<-length(s)
  freq_a<-count_A/m
  freq_t<-count_T/m
  freq_g<-count_G/m
  freq_c<-count_C/m
  freq<-data.frame(a_n=a_n,a=freq_a,t=freq_t,g=freq_g,c=freq_c)
  print(freq)
}
x<-read_genbank("AB115403.gb")
x

#function2
genbank2fasta<-function(file){
  seq<-readLines(file)
  a_n<-paste0(">",unlist(strsplit(seq[1],split = "\\s+"))[2])
  st<-grep("ORIGIN",seq)
  ed<-grep("^//",seq)
  gb_seq<-gsub(" ","",gsub("\\d","",paste(seq[(st+1):(ed-1)],collapse = "")))
  res<-paste0(a_n,"\n",gb_seq)
  return(res)
}
y<-genbank2fasta("AB115403.gb")
y
writeLines(y)
cat(y)
write.table(y,
            file = "AB115403.fas",
            row.names = F,
            col.names = F,
            quote = F)

#if we have 2 or more sequences
x<-genbank2fasta("AJ534526.gb")
y<-genbank2fasta("AB115403.gb")
z<-paste0(x,"\n",y)
writeLines(z)
cat(z)
write.table(z,
            file = "test.fas",
            row.names = F,
            col.names = F,
            quote = F)

###another methods
#BiocManager::install("genbankr")
#library(genbankr)
#gb1<-readGenBank("AB115403.gb")
library(ape)
gb2<-read.GenBank("AB115403",as.character = T)

###attention
read.genbank <- function(file) {
  x <- readLines(file)
  acc <- sub("\\w+\\s+(\\w+)$", "\\1", x[grep("^ACCESSION", x)])
  src <- sub("SOURCE\\s+",  "", x[grep("^SOURCE",x)])
  header <- paste0(src, "(", acc, ")")
  i <- grep("ORIGIN", x)
  ss <- x[(i+1):length(x)]
  ss <- ss[1:(grep("//", ss) -1)]
  ss <- gsub("\\s+\\d+", "", ss)
  ss <- gsub("\\s+", "", ss)
  ss <- paste0(ss, collapse = "")
  f <- tempfile(fileext = ".fasta")
  cat(">",  file = f,append = TRUE)
  cat(header,  file = f,append = TRUE, sep = "\n")
  cat(ss,  file = f,append = TRUE, sep = "\n")
  read.fasta(f)
}
write.fasta <- function(x, file) {
  ape::write.FASTA(x, file)
}
x<-read.genbank("AB115403.gb")
f<-tempfile(fileext = ".fa")
write.fasta(x,f)
readLines(f)
read.fasta(f)