#################################################################
#model start for pfos (adme), PC,Vol, blood flow


library(pksensi)
library(sensitivity)

library(mcmcr)
# theme_set(theme_light())

setwd("~/MCSim/mod")
library("ggformula")

Modeling_dir = "NFT_PBPK2"
directory_i = "~/MCSim/mod/NFT_PBPK2/Manuscript_figure"


makemcsim <- function(model,
                      deSolve = F,
                      dir = Modeling_dir) {
  exe_file <- paste0("mcsim_", model, " ")
  
  if (file.exists(Modeling_dir) == F) {
    dir.create(Modeling_dir)
  }
  
  if (file.exists(model)) {
    # Move model file from working directory to Modeling folder
    invisible(file.copy(
      from = paste0(getwd(), "/", model),
      to = paste0(getwd(), "/",Modeling_dir,"/", model)
    ))
    invisible(file.remove(model))
    message(paste0("* The '", model, "' has been moved to Modeling folder."))
  }
  
  if (deSolve == T) {
    system(paste("./MCSim/mod -R ", dir, "/", model, " ", model, ".c", sep = ""))
    system (paste0("R CMD SHLIB ", model, ".c")) # create *.dll files
    dyn.load(paste(model, .Platform$dynlib.ext, sep = "")) # load *.dll
    source(paste0(model, "_inits.R"))
  } else {
    system(paste("../mod/mod ", dir, "/", model, " ", model, ".c", sep = ""))
    system(
      paste(
        "gcc -O3 -I.. -I../sim -o mcsim_",
        model,
        " ",
        model,
        ".c ../sim/*.c -lm ",
        sep = ""
      )
    )
    invisible(file.remove(paste0(model, ".c")))
    if (file.exists(exe_file))
      message(paste0("* Created executable program '", exe_file, "'."))
  }
}

#source("C:/Users/alexe/OneDrive - Universiteit Leiden/MCSim/mod/Mcsim_source_function.R") ###read the function file directly

####################################

model <- "Rabbit_V5asensitivity.model.R"
makemcsim(model = model, dir = Modeling_dir)  #change the compiled file name replacing mcsim_filename to mcsim.filename fo the file 

Parmeters = read.csv(paste0(Modeling_dir,"/","Posterior_chains/","RabbitV5mixed_NFTPBPK.csv"))


# create a vector of the values from the CSV file
params_vector <- c(Parmeters$Value)

# name the vector using the names from the CSV file
names(params_vector) <- Parmeters$Parms

setnames = Parmeters$Parms
# assign the named vector to the variable 'parms'
parms <- params_vector


parm2 = parms[1:length(parms)]*(1 - 0.5)
parm3 = parms[1:length(parms)]*(1 +1.5)

dist <- rep("Uniform", length(parms))
q <- rep("qunif", length(parms))

q.arg <- list()

# loop through the elements of parm2 and parm3
for (i in 1:length(parms)) {
  # use paste to create the list of lists
  q.arg[[i]] <- list(parm2[[i]], parm3[[i]])
}



#generate parameter matrix
set.seed(1234)

vars <- c("cplasma","ckidney","cliver", "cgut", "crestbody", "ctubules","Aurine")

#vars <- c("cplasma")
times <- seq(from = 1, to = 24, by = 1)

x <- rfast99(params = setnames, n = 200, q = q, q.arg = q.arg, replicate = 1)


#conditions = c("oraldose1=PerDose(20,24,0,0.001)", "oraldose= 20");

conditions = c("IVdosing=15", "oraldose= 0");

#go change the file name by remvoing underscore in mod to dot.
y_SUL <- solve_mcsim(x, mName = model,
                     params = setnames,
                     time = times,
                     vars = vars,
                     setpoint.name = setnames,
                     condition = conditions)

tiff(paste0(directory_i,"/", "FigS12_15mg_kg_rabbitIV",".png"), units="in", width=10, height=8, res=700, compression = 'lzw')

par(mfrow = c(4,5),
    oma = c(5,4,1,1) + 0.1,
    mar = c(1,1,1,1) + 1)
P = heat_check(y_SUL, times = c(5, 10, 15, 20),order = "total order", show.all = T)
print(P)
dev.off()

Y <- solve_mcsim(x, mName = model,
                     params = setnames, 
                     time = times, 
                     vars = vars,
                     condition = conditions, 
                     rtol = 1e-7, atol = 1e-9)

plot(Y)
# 
# heat_check(Y, order = "total order", show.all = T)
# 
# heat_check(Y, index = "CI", order = "total order", show.all = T)

final_tSI = cbind(y_SUL$tSI[,1:length(parms),1])

final_mSI = cbind(y_SUL$mSI[,1:length(parms),1])


final_iSI = cbind(y_SUL$iSI[,1:length(parms),1])

# colnames(y_Sul32R$tSI) = c("","kd","Kmu" ,"nm", 
#                            "Vmax2","Kmm","nu","r1.Sul",'d1.Sul',"r2.Sul32",'d2.Sul32')


figdirec = "~/Modeling_dir"
tiff(paste0(figdirec,"/", "Sensitivity_PBPKNFT",".png"), units="in", width=10, height=8, res=700, compression = 'lzw')
# pdf(paste0(figdirec,"/", "Sensitivity1",".pdf"),onefile=FALSE)
par(mfrow = c(4,8))
par(mfrow = c(4,5),
    oma = c(5,4,1,1) + 0.1,
    mar = c(1,1,1,1) + 1)
for (i in 1:length(parms)){
  plot(times, final_tSI[,i],ylim = c(0, 0.9), bty = "n", yaxt = "n",type = "l", col = "black",main =colnames(final_tSI)[i],ylab = "SI",lwd = 1.5,cex.main=1.0,cex.lab=1.0, cex.axis=1.0,xlab="Time[h]")
  #lines(times, final_mSI_R48_32[,i], col = "red",ylim = c(0, 1),yaxt = "n",lwd = 1.5)
  lines(times, final_iSI[,i], col = "red",ylim = c(0, 0.9),yaxt = "n",lty = 3,lwd = 1.5)
  abline(h = 0.05, col = "blue", lty = 2,lwd = 1.2)
  axis(2, at = seq(0, 1, .2))
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom',legend = c("Total order", "first order","cut-off"), col = c("black","red", "blue"), lty=c(1,1,2),lwd = 2, xpd = TRUE, horiz = TRUE, cex = 1, seg.len=2, bty = 'n')

# legend('bottom',legend = c("Total order", "first order", "mixed order", "cut-off"), col = c("black", "grey", "red", "blue"), lty=c(1,3,1,2),lwd = 2, xpd = TRUE, horiz = TRUE, cex = 1, seg.len=2, bty = 'n')
# # xpd = TRUE makes the legend plot to the figure
dev.off()





par(mfrow = c(5,5))
for (i in 1:length(parms)){
  plot(times, final_tSI[,i],ylim = c(0, 1), bty = "n", yaxt = "n",type = "l", col = "black",main =colnames(final_tSI)[i],ylab = "SI",lwd = 2)
  lines(times, final_mSI[,i], col = "red",ylim = c(0, 1),yaxt = "n",lwd = 2)
  abline(h = 0.1, col = "black", lty = 2)
  axis(2, at = seq(0, 1, .2))
}
