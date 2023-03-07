library("rlang")
library("sodium")
library("devtools")
library("shinymanager")
library("DT")
library("dplyr")
library("fontawesome")
library("shinydashboard")
library("shinycssloaders")
library("igraph")
library("rmarkdown")
library("markdown")
library("ggplot2")
library("visNetwork")
library("heatmaply")
library("plotly")
library("brew")
library("RColorBrewer")
library("curl")
library("wesanderson")
library("shinyauthr")

rm(list=ls())
gc()
# options(repos = BiocManager::repositories())
#.libPaths("/home/mbotos/R/x86_64-pc-linux-gnu-library/4.2/")


inactivity <- "function idleTimer() {
var t = setTimeout(logout, 120000);
window.onmousemove = resetTimer; // catches mouse movements
window.onmousedown = resetTimer; // catches mouse movements
window.onclick = resetTimer;     // catches mouse clicks
window.onscroll = resetTimer;    // catches scrolling
window.onkeypress = resetTimer;  //catches keyboard actions

function logout() {
window.close();  //close the window
}

function resetTimer() {
clearTimeout(t);
t = setTimeout(logout, 120000);  // time is in milliseconds (1000 is 1 second)
}
}
idleTimer();"


# data.frame with credentials info
credentials <- data.frame(
  user = c("marius","nadia","panos","prateek", ""),
  password = c("marius","nadia","panos","prateek", ""),
  start = c("2022-04-01","2022-04-01","2022-04-01","2022-04-01", "2022-04-01"),
  expire = c(NA,NA,NA,NA, NA),
  admin = c(TRUE,FALSE,FALSE,TRUE, TRUE),
  comment = "Simple and secure authentification mechanism
  for single ‘Shiny’ applications.",
  stringsAsFactors = FALSE
)




#################################
#                               #
#   Data Load                   #
#                               #     
#################################
#     Ablated                   #
#################################
##:%s/^\(.*\)$/"\1",/
#cat NN | tr "\n" " "
#excel to paste and make the final command
#setwd("/media/marius/Samsung_T5/PhD/Projects/prsa/Outputs/r/Cytoscape/Cleaned_Annotation_Ensembl_v4/ResultsCytoscape/PossibleShiny/RShinyCode/ShinnyNetworksGOBP")


myCol <- c("#CC79A7","#0072B2","#E69F00","#999999")
myCol_3 <- c("#CC79A7","#0072B2","#009E73","#999999")


abl <- read_graph(file = "./data/cleaned_unique_go_bp_v2__Ablated_vs_Sham.graphml",format = "graphml")

npp <- c("GO:1901293", "GO:0072522", "GO:0006163", "GO:0022900", "GO:0042775", "GO:0006123", "GO:1903715", "GO:0046040", "GO:0034033", "GO:0033866", "GO:0015986", "GO:0006085", "GO:0044272", "GO:0072521", "GO:0043467", "GO:0072350", "GO:0006188", "GO:0043648", "GO:0071616", "GO:0009205", "GO:0042773", "GO:0009201", "GO:0006084", "GO:0009142", "GO:0019646", "GO:0034032", "GO:0006753", "GO:0019693", "GO:0006790", "GO:0033865", "GO:0009145", "GO:0046390", "GO:0009127", "GO:0009199", "GO:0006164", "GO:0009152", "GO:0035383", "GO:0033875", "GO:0009060", "GO:0019359", "GO:0006637", "GO:0009260", "GO:0015980", "GO:0009124", "GO:1902600", "GO:0006167", "GO:0006099", "GO:0009144", "GO:0009259", "GO:0006107", "GO:0006086", "GO:0015985", "GO:0009165", "GO:0045333", "GO:0034030", "GO:0046034", "GO:0035384", "GO:0042776", "GO:0009156", "GO:0006120", "GO:0006754", "GO:0022904", "GO:0046033", "GO:0009206", "GO:0006090", "GO:0009117", "GO:0006119", "GO:0009141", "GO:0043457", "GO:0009168")
cvt <- c("GO:0055010", "GO:0060415", "GO:0003228", "GO:0060343", "GO:0003206", "GO:0048644", "GO:0003208", "GO:0003231", "GO:0003209", "GO:0003230", "GO:0055008", "GO:0003205", "GO:0061383", "GO:0003214", "GO:0055015", "GO:0048845", "GO:0061384", "GO:0001570", "GO:0055012", "GO:0003229", "GO:0003007", "GO:0060841", "GO:0060347", "GO:0003222")
ash <- c("GO:0046165", "GO:0042632", "GO:0008209", "GO:0042572", "GO:0006081", "GO:0042446", "GO:0016125", "GO:0034754", "GO:0008210", "GO:0032365", "GO:0032366", "GO:0008203", "GO:0034308", "GO:0055088", "GO:1902652", "GO:0032367", "GO:0042445", "GO:0055092", "GO:0006066")
gpc <- c("GO:0010675", "GO:0010676", "GO:0032881", "GO:0001890", "GO:0010906", "GO:0045913", "GO:0009250", "GO:0033692", "GO:0010907", "GO:0005996", "GO:0005978", "GO:0070875", "GO:0034637", "GO:0043255", "GO:0016051", "GO:0019318", "GO:0000271", "GO:0006109", "GO:0006006")
erd <- c("GO:0045685", "GO:0014014", "GO:0045604", "GO:0042491", "GO:0045631", "GO:0014013", "GO:0010001", "GO:0061448", "GO:0042490", "GO:2000980", "GO:0045686", "GO:0030856", "GO:0045607", "GO:0048710", "GO:0050768")
oaf <- c("GO:0044242", "GO:0046322", "GO:0019395", "GO:0006631", "GO:0033539", "GO:0006575", "GO:0034440", "GO:0009437", "GO:0006635", "GO:0046459", "GO:0072329", "GO:0016042", "GO:0051791", "GO:0006577", "GO:0009062")
nma <- c("GO:0014733", "GO:0043500", "GO:0014741", "GO:0014888", "GO:0043501", "GO:0050879", "GO:0003009", "GO:0010614", "GO:0014883", "GO:0050881", "GO:0050885", "GO:0050905", "GO:0043502")
vsa <- c("GO:0045823", "GO:0035904", "GO:0097084", "GO:2001212", "GO:0060840", "GO:1903524", "GO:0035886", "GO:0060947", "GO:0010460", "GO:0003184", "GO:0048844", "GO:0002026", "GO:0009798")
dms <- c("GO:0072673", "GO:0050773", "GO:0022604", "GO:0010770", "GO:0008360", "GO:0060997", "GO:0050775", "GO:0048813", "GO:0048814", "GO:0016358", "GO:0060996", "GO:0010769")
mpm <- c("GO:0097345", "GO:1905710", "GO:0006839", "GO:1902110", "GO:1902686", "GO:0090559", "GO:0035794", "GO:0008637", "GO:0007006", "GO:0046902", "GO:1902108", "GO:0010821")
ota <- c("GO:0098656", "GO:0003333", "GO:0015816", "GO:0015804", "GO:0015807", "GO:1905039", "GO:0015711", "GO:0015849", "GO:0043090", "GO:0006865", "GO:0032328", "GO:1903825")
#
#vma <- vertex_attr(abl)$name[vertex_attr(abl)$`__glayCluster` == 1]
vma <- c("GO:0072384", "GO:0008089", "GO:0051650", "GO:0008090", "GO:0047496", "GO:0008088", "GO:0098930", "GO:0051648", "GO:0010970", "GO:0006900", "GO:0099518")
its <- c("GO:0050729", "GO:0031620", "GO:0044070", "GO:0031650", "GO:0002673", "GO:0002675", "GO:0001660", "GO:0002526", "GO:0009636", "GO:0031649")
smm <- c("GO:0006582", "GO:0019748", "GO:0048021", "GO:0043455", "GO:0042438", "GO:1900376", "GO:0042440", "GO:0044550", "GO:0046148", "GO:0046189")
mml <- c("GO:0002687", "GO:0032103", "GO:1905521", "GO:0050920", "GO:0050900", "GO:1905517", "GO:0060326", "GO:0097529", "GO:0002685")
pca <- c("GO:0006576", "GO:0097164", "GO:0042417", "GO:0033240", "GO:0042428", "GO:0042069", "GO:0009712", "GO:1901160", "GO:0006584")
scg <- c("GO:0006672", "GO:0006677", "GO:0046479", "GO:0001573", "GO:0019377", "GO:0030149", "GO:0046466", "GO:0006664", "GO:1903509")
dmv <- c("GO:0007189", "GO:0097746", "GO:0035296", "GO:0046541", "GO:0019229", "GO:0042311", "GO:0042310", "GO:0035150")
eem <- c("GO:0010595", "GO:0043542", "GO:0035767", "GO:0010634", "GO:0010594", "GO:2001026", "GO:0043536", "GO:0010632")
kca <- c("GO:0007254", "GO:0051403", "GO:0043406", "GO:0071902", "GO:0043507", "GO:0043405", "GO:0043506", "GO:0031098")
mei <- c("GO:0031623", "GO:0045807", "GO:0030100", "GO:0048259", "GO:0002090", "GO:0043112", "GO:0032801", "GO:0006898")
ppd <- c("GO:0016050", "GO:0032438", "GO:0033059", "GO:0048070", "GO:0048066", "GO:0120305", "GO:0050931", "GO:0048753")
pga <- c("GO:0006024", "GO:1903510", "GO:0015012", "GO:0030203", "GO:0006023", "GO:0030202", "GO:0006022", "GO:0030166")
sag <- c("GO:0002040", "GO:0001569", "GO:0040036", "GO:1903587", "GO:0044344", "GO:0002043", "GO:0008543", "GO:0071774")
tvl <- c("GO:0006623", "GO:0007041", "GO:0061462", "GO:0006605", "GO:0071806", "GO:0006622", "GO:0072665", "GO:0072594")
dsc <- c("GO:0051937", "GO:0015872", "GO:0050432", "GO:0014046", "GO:0050433", "GO:0014059", "GO:0015844")
kqu <- c("GO:1901663", "GO:1901661", "GO:0006743", "GO:0042181", "GO:0042180", "GO:0042182", "GO:0006744")
pnc <- c("GO:0021953", "GO:0021954", "GO:0021872", "GO:0021761", "GO:0021879", "GO:0021860", "GO:0021859")
cpc <- c("GO:0060920", "GO:0055006", "GO:0035051", "GO:0055013", "GO:0055007", "GO:0060926")
ccr <- c("GO:0021801", "GO:0022030", "GO:0008347", "GO:0021799", "GO:0042063", "GO:0021819")
ela <- c("GO:0035113", "GO:0035107", "GO:0030326", "GO:0048736", "GO:0060173", "GO:0035108")
ese <- c("GO:0042058","GO:0007173","GO:0045742","GO:1901186","GO:1901184","GO:0038127")
ehf <- c("GO:0098773","GO:0022405","GO:0022404","GO:0008544","GO:0048730","GO:0031069")
gfs <- c("GO:0007548","GO:0045137","GO:0048806","GO:0035112","GO:0030539","GO:0046660")
pgp <- c("GO:0046474","GO:0046486","GO:0008654","GO:0006661","GO:0034638","GO:0071071")
ros <- c("GO:1903426","GO:1903409","GO:1903427","GO:0072593","GO:2000377","GO:2000378")
srg <- c("GO:0050810","GO:0034698","GO:0045939","GO:0071371","GO:0019218","GO:0006694")
sua <- c("GO:0009225","GO:0006047","GO:0006044","GO:0006040","GO:0006054","GO:0009226")
sdb <- c("GO:0009855","GO:0035050","GO:0060972","GO:0007368","GO:0060973","GO:0009799")
sde <- c("GO:0016081","GO:0048278","GO:0022406","GO:0006904","GO:0140029","GO:0140056")
tic <- c("GO:0001659","GO:0120162","GO:0002024","GO:0106106","GO:0120161","GO:1990845")
cim <- c("GO:1990542","GO:0090279","GO:0006851","GO:0070509","GO:0036444")
cel <- c("GO:1905953","GO:1905952","GO:0033344","GO:0010874","GO:0032369")
mtg <- c("GO:0034219","GO:0010827","GO:0010828","GO:0015749","GO:0008643")
mcc <- c("GO:0060048","GO:0006941","GO:0006942","GO:0055117","GO:1903779")
ncs <- c("GO:0048863","GO:0014032","GO:0014033","GO:0001755","GO:0048864")
put <- c("GO:0051438","GO:0031396","GO:0051444","GO:2000434","GO:1903320")
rrl <- c("GO:0001974","GO:0030324","GO:0030323","GO:0060541","GO:0001886")
rlh <- c("GO:0071453","GO:0036293","GO:0071456","GO:0070482","GO:0001666")
scp <- c("GO:0014841","GO:0014842","GO:0014857","GO:0014855","GO:0014856")
sip <- c("GO:0006814","GO:0030007","GO:0002028","GO:0055075","GO:0030004")
abw <- c("GO:0090659","GO:0008344","GO:0030534","GO:0007628")
bra <- c("GO:0045453","GO:0060249","GO:0001894","GO:0046849")
dca <- c("GO:0017004","GO:0032981","GO:0033108","GO:0010257")
get <- c("GO:0006390","GO:0140053","GO:0032543","GO:0062125")
ias <- c("GO:2001243","GO:2001233","GO:2001234","GO:1902236")
lpp <- c("GO:0072659","GO:1905475","GO:1904375","GO:1990778")
mas <- c("GO:0010665","GO:0010662","GO:0010661","GO:0010658")
mfr <- c("GO:0014902","GO:0010830","GO:1901739","GO:0007520")
msc <- c("GO:0045214","GO:0030239","GO:0010927","GO:0055002")
pnr <- c("GO:0048755","GO:0070572","GO:0050772","GO:0048680")
rsc <- c("GO:0070296","GO:1903514","GO:0014808","GO:0010880")
tgg <- c("GO:0007215","GO:0051966","GO:0035249","GO:0007216")
tcc <- c("GO:2001259","GO:0032411","GO:0032414","GO:1901021")
tpp <- c("GO:0050731","GO:0050730","GO:0018108","GO:0018212")
bbk <- c("GO:0070857","GO:0006699","GO:0010565")
cmp <- c("GO:0042407","GO:0007007","GO:0034982")
cgm <- c("GO:0044839","GO:0000086","GO:1902749")
cte <- c("GO:0043154","GO:2000117","GO:0043281")
enc <- c("GO:0002262","GO:0048872","GO:0034101")
gra <- c("GO:0001975","GO:0014075","GO:0006177")
iam <- c("GO:0014072","GO:0043278","GO:0043279")
lljs <- c("GO:0035418","GO:1902414","GO:1902473")
mgl <- c("GO:0006486","GO:0043413","GO:0006487")
mog <- c("GO:0040015","GO:0035264","GO:0040014")
nf6 <- c("GO:0061615","GO:0016052","GO:0019674")
nlc <- c("GO:0045833","GO:0050994","GO:0050995")
nna <- c("GO:0051402","GO:0043524","GO:0043523")
nkt <- c("GO:0051092","GO:0051090","GO:0043433")
nsp <- c("GO:0008593","GO:0045747","GO:0007219")
ppt <- c("GO:0018107","GO:0018105","GO:0018210")
paj <- c("GO:0051965","GO:1901890","GO:0099560")
pvp <- c("GO:0043116","GO:0043117","GO:0043114")
psd <- c("GO:0097106","GO:0099558","GO:0099084")
sgs <- c("GO:0033500","GO:0042593","GO:0061179")
shi <- c("GO:0071383","GO:0030518","GO:0033143")
vef <- c("GO:0030949","GO:0048010","GO:0030947")
api <- c("GO:1902742","GO:0060561")
bca <- c("GO:0009083","GO:0009081")
cmh <- c("GO:0030902","GO:0021535")
ci7 <- c("GO:0098760","GO:0098761")
fao <- c("GO:0046320","GO:0045834")
hde <- c("GO:0061008","GO:0003157")
irp <- c("GO:0043567","GO:0048009")
ler <- c("GO:0072599","GO:0070972")
mmp <- c("GO:0010918","GO:0045838")
mab <- c("GO:0072330","GO:0006633")
nrc <- c("GO:0042754","GO:0019915")
odt <- c("GO:0042476","GO:0042475")
pla <- c("GO:0051349","GO:0051339")
poa <- c("GO:0051341","GO:0051353")
pmo <- c("GO:0043113","GO:0001941")
pks <- c("GO:0051896","GO:0043491")
rgc <- c("GO:0007178","GO:0060019")
ross <- c("GO:0034614","GO:0000302")
rcr <- c("GO:0007623","GO:0042753")
rfp <- c("GO:0048144","GO:0048145")
rla <- c("GO:0007159","GO:1903037")
rmp <- c("GO:2000291","GO:0051450")
rcm <- c("GO:0055119","GO:0090075")
rge <- c("GO:0003016","GO:0007585")
rtc <- c("GO:0009266","GO:0009409")
sss <- c("GO:0007379","GO:0035282")
spp <- c("GO:0016486","GO:0140448")
vda <- c("GO:1904018","GO:0045766")


vertex_attr(abl)$AA <- case_when(vertex_attr(abl)$name %in%	npp 	~ "nucleoside purine process",
                                 vertex_attr(abl)$name %in%	cvt 	~ "cardiac ventricle trabecula",
                                 vertex_attr(abl)$name %in%	ash 	~ "alcohol sterol homeostasis",
                                 vertex_attr(abl)$name %in%	gpc 	~ "glycogen polysaccharide carbohydrate",
                                 vertex_attr(abl)$name %in%	erd 	~ "ear receptor differentiation",
                                 vertex_attr(abl)$name %in%	oaf 	~ "oxidation acid fatty",
                                 vertex_attr(abl)$name %in%	nma 	~ "neuromuscular movement adaptation",
                                 vertex_attr(abl)$name %in%	vsa 	~ "vascular smooth artery",
                                 vertex_attr(abl)$name %in%	dms 	~ "dendrite morphogenesis spine",
                                 vertex_attr(abl)$name %in%	mpm 	~ "membrane permeability mitochondrial",
                                 vertex_attr(abl)$name %in%	ota 	~ "organic transport amino",
                                 vertex_attr(abl)$name %in%	vma 	~ "vesicle microtubule axonal",
                                 vertex_attr(abl)$name %in%	its 	~ "inflammatory toxic substance",
                                 vertex_attr(abl)$name %in%	smm 	~ "secondary metabolite melanin",
                                 vertex_attr(abl)$name %in%	mml 	~ "macrophage migration leukocyte",
                                 vertex_attr(abl)$name %in%	pca 	~ "positive cellular amine",
                                 vertex_attr(abl)$name %in%	scg 	~ "sphingolipid catabolic glycolipid",
                                 vertex_attr(abl)$name %in%	dmv 	~ "diameter maintenance vasoconstriction",
                                 vertex_attr(abl)$name %in%	eem 	~ "endothelial epithelial migration",
                                 vertex_attr(abl)$name %in%	kca 	~ "kinase cascade activated",
                                 vertex_attr(abl)$name %in%	mei 	~ "mediated endocytosis internalization",
                                 vertex_attr(abl)$name %in%	ppd 	~ "pigmentation pigment developmental",
                                 vertex_attr(abl)$name %in%	pga 	~ "proteoglycan glycosaminoglycan aminoglycan",
                                 vertex_attr(abl)$name %in%	sag 	~ "sprouting angiogenesis fibroblast",
                                 vertex_attr(abl)$name %in%	tvl 	~ "targeting vacuole lysosome",
                                 vertex_attr(abl)$name %in%	dsc 	~ "dopamine secretion catecholamine",
                                 vertex_attr(abl)$name %in%	kqu 	~ "ketone quinone ubiquinone",
                                 vertex_attr(abl)$name %in%	pnc 	~ "pyramidal neuron central",
                                 vertex_attr(abl)$name %in%	cpc 	~ "cardiac pacemaker cardiocyte",
                                 vertex_attr(abl)$name %in%	ccr 	~ "cerebral cortex radially",
                                 vertex_attr(abl)$name %in%	ela 	~ "embryonic limb appendage",
                                 vertex_attr(abl)$name %in%	ese 	~ "epidermal signaling erbb",
                                 vertex_attr(abl)$name %in%	ehf 	~ "epidermis hair follicle",
                                 vertex_attr(abl)$name %in%	gfs 	~ "genitalia female sex",
                                 vertex_attr(abl)$name %in%	pgp 	~ "phospholipid glycerolipid phosphatidylinositol",
                                 vertex_attr(abl)$name %in%	ros 	~ "reactive oxygen species",
                                 vertex_attr(abl)$name %in%	srg 	~ "steroid response gonadotropin",
                                 vertex_attr(abl)$name %in%	sua 	~ "sugar udp acetylglucosamine",
                                 vertex_attr(abl)$name %in%	sdb 	~ "symmetry determination bilateral",
                                 vertex_attr(abl)$name %in%	sde 	~ "synaptic docking exocytosis",
                                 vertex_attr(abl)$name %in%	tic 	~ "thermogenesis induced cold",
                                 vertex_attr(abl)$name %in%	cim 	~ "calcium import mitochondrion",
                                 vertex_attr(abl)$name %in%	cel 	~ "cholesterol efflux lipid",
                                 vertex_attr(abl)$name %in%	mtg 	~ "monosaccharide transmembrane glucose",
                                 vertex_attr(abl)$name %in%	mcc 	~ "muscle contraction conduction",
                                 vertex_attr(abl)$name %in%	ncs 	~ "neural crest stem",
                                 vertex_attr(abl)$name %in%	put 	~ "protein ubiquitin transferase",
                                 vertex_attr(abl)$name %in%	rrl 	~ "respiratory remodeling lung",
                                 vertex_attr(abl)$name %in%	rlh 	~ "response levels hypoxia",
                                 vertex_attr(abl)$name %in%	scp 	~ "satellite skeletal proliferation",
                                 vertex_attr(abl)$name %in%	sip 	~ "sodium ion potassium",
                                 vertex_attr(abl)$name %in%	abw 	~ "adult behavior walking",
                                 vertex_attr(abl)$name %in%	bra 	~ "bone resorption anatomical",
                                 vertex_attr(abl)$name %in%	dca 	~ "dehydrogenase complex assembly",
                                 vertex_attr(abl)$name %in%	get 	~ "gene expression translation",
                                 vertex_attr(abl)$name %in%	ias 	~ "intrinsic apoptotic signaling",
                                 vertex_attr(abl)$name %in%	lpp 	~ "localization plasma periphery",
                                 vertex_attr(abl)$name %in%	mas 	~ "muscle apoptotic striated",
                                 vertex_attr(abl)$name %in%	mfr 	~ "myoblast fusion regulation",
                                 vertex_attr(abl)$name %in%	msc 	~ "myofibril sarcomere component",
                                 vertex_attr(abl)$name %in%	pnr 	~ "positive neuron regeneration",
                                 vertex_attr(abl)$name %in%	rsc 	~ "release sequestered cytosol",
                                 vertex_attr(abl)$name %in%	tgg 	~ "transmission glutamatergic glutamate",
                                 vertex_attr(abl)$name %in%	tcc 	~ "transporter cation channel",
                                 vertex_attr(abl)$name %in%	tpp 	~ "tyrosine peptidyl phosphorylation",
                                 vertex_attr(abl)$name %in%	bbk 	~ "bile biosynthetic ketone",
                                 vertex_attr(abl)$name %in%	cmp 	~ "cristae mitochondrial processing",
                                 vertex_attr(abl)$name %in%	cgm 	~ "cycle g2 mitotic",
                                 vertex_attr(abl)$name %in%	cte 	~ "cysteine type endopeptidase",
                                 vertex_attr(abl)$name %in%	enc 	~ "erythrocyte number cells",
                                 vertex_attr(abl)$name %in%	gra 	~ "gmp response amphetamine",
                                 vertex_attr(abl)$name %in%	iam 	~ "isoquinoline alkaloid morphine",
                                 vertex_attr(abl)$name %in%	lljs 	~ "localization junction synapse",
                                 vertex_attr(abl)$name %in%	mgl 	~ "macromolecule glycosylation linked",
                                 vertex_attr(abl)$name %in%	mog 	~ "multicellular organism growth",
                                 vertex_attr(abl)$name %in%	nf6 	~ "nad fructose 6",
                                 vertex_attr(abl)$name %in%	nlc 	~ "negative lipid catabolic",
                                 vertex_attr(abl)$name %in%	nna 	~ "negative neuron apoptotic",
                                 vertex_attr(abl)$name %in%	nkt 	~ "nf kappab transcription",
                                 vertex_attr(abl)$name %in%	nsp 	~ "notch signaling pathway",
                                 vertex_attr(abl)$name %in%	ppt 	~ "peptidyl phosphorylation threonine",
                                 vertex_attr(abl)$name %in%	paj 	~ "positive assembly junction",
                                 vertex_attr(abl)$name %in%	pvp 	~ "positive vascular permeability",
                                 vertex_attr(abl)$name %in%	psd 	~ "postsynaptic specialization density",
                                 vertex_attr(abl)$name %in%	sgs 	~ "secretion glucose stimulus",
                                 vertex_attr(abl)$name %in%	shi 	~ "steroid hormone intracellular",
                                 vertex_attr(abl)$name %in%	vef 	~ "vascular endothelial factor",
                                 vertex_attr(abl)$name %in%	api 	~ "apoptotic process involved",
                                 vertex_attr(abl)$name %in%	bca 	~ "branched chain acid",
                                 vertex_attr(abl)$name %in%	cmh 	~ "cell migration hindbrain",
                                 vertex_attr(abl)$name %in%	ci7 	~ "cellular interleukin 7",
                                 vertex_attr(abl)$name %in%	fao 	~ "fatty acid oxidation",
                                 vertex_attr(abl)$name %in%	hde 	~ "hepaticobiliary development endocardium",
                                 vertex_attr(abl)$name %in%	irp 	~ "insulin receptor pathway",
                                 vertex_attr(abl)$name %in%	ler 	~ "localization endoplasmic reticulum",
                                 vertex_attr(abl)$name %in%	mmp 	~ "mitochondrial membrane potential",
                                 vertex_attr(abl)$name %in%	mab 	~ "monocarboxylic acid biosynthetic",
                                 vertex_attr(abl)$name %in%	nrc 	~ "negative regulation circadian",
                                 vertex_attr(abl)$name %in%	odt 	~ "odontogenesis dentin tooth",
                                 vertex_attr(abl)$name %in%	pla 	~ "positive lyase activity",
                                 vertex_attr(abl)$name %in%	poa 	~ "positive oxidoreductase activity",
                                 vertex_attr(abl)$name %in%	pmo 	~ "postsynaptic membrane organization",
                                 vertex_attr(abl)$name %in%	pks 	~ "protein kinase signaling",
                                 vertex_attr(abl)$name %in%	rgc 	~ "radial glial cell",
                                 vertex_attr(abl)$name %in%	ross 	~ "reactive oxygen species",
                                 vertex_attr(abl)$name %in%	rcr 	~ "regulation circadian rhythm",
                                 vertex_attr(abl)$name %in%	rfp 	~ "regulation fibroblast proliferation",
                                 vertex_attr(abl)$name %in%	rla 	~ "regulation leukocyte adhesion",
                                 vertex_attr(abl)$name %in%	rmp 	~ "regulation myoblast proliferation",
                                 vertex_attr(abl)$name %in%	rcm 	~ "relaxation cardiac muscle",
                                 vertex_attr(abl)$name %in%	rge 	~ "respiratory gaseous exchange",
                                 vertex_attr(abl)$name %in%	rtc 	~ "response temperature cold",
                                 vertex_attr(abl)$name %in%	sss 	~ "segment specification segmentation",
                                 vertex_attr(abl)$name %in%	spp 	~ "signaling precursor processing",
                                 vertex_attr(abl)$name %in%	vda	~ "vasculature development angiogenesis")


long_list_to_fill_abl <- lapply(seq_along(unique(vertex_attr(abl)$AA)),function(i) i)

AA_Clusters <- lapply(seq_along(vertex_attr(abl)$AA),function(n) {
  if_else(is.na(vertex_attr(abl)$AA[[n]]) == TRUE,
          long_list_to_fill_abl[n] <- vertex_attr(abl)$`EnrichmentMap::GS_DESCR`[[n]],
          long_list_to_fill_abl[n] <- vertex_attr(abl)$AA[[n]])
})

vertex_attr(abl)$Clusters <- unlist(AA_Clusters)

#################################
#  Resection                    #
#################################
amp <- read_graph(file = "./data/cleaned_unique_go_bp_v2__Resection_vs_Sham.graphml",
                  format = "graphml")

pim <- c("GO:0001819", "GO:0002263", "GO:0002274", "GO:0002275", "GO:0002279", "GO:0002366", "GO:0002420", "GO:0002444", "GO:0002448", "GO:0002697", "GO:0002699", "GO:0002700", "GO:0002702", "GO:0002703", "GO:0002705", "GO:0002717", "GO:0002718", "GO:0002858", "GO:0002886", "GO:0031343", "GO:0032418", "GO:0033003", "GO:0033005", "GO:0033006", "GO:0033008", "GO:0043299", "GO:0043300", "GO:0043302", "GO:0043303", "GO:0043304", "GO:0043306", "GO:0045576", "GO:0045954", "GO:0051656")
brd <- c("GO:0000018", "GO:0000724", "GO:0000725", "GO:0000727", "GO:0000731", "GO:0006260", "GO:0006261", "GO:0006268", "GO:0006270", "GO:0006275", "GO:0006281", "GO:0006284", "GO:0006302", "GO:0006310", "GO:0006336", "GO:0010216", "GO:0032392", "GO:0032508", "GO:0033260", "GO:0034724", "GO:0044786", "GO:0045740", "GO:0071103", "GO:0090329", "GO:1900262", "GO:1900264", "GO:2000042")

pcl <- c("GO:0001783", "GO:0002696", "GO:0002902", "GO:0002903", "GO:0022407", "GO:0022409", "GO:0030217", "GO:0030888", "GO:0030890", "GO:0032943", "GO:0032944", "GO:0032946", "GO:0033077", "GO:0042100", "GO:0043368", "GO:0045058", "GO:0045785", "GO:0046631", "GO:0046651", "GO:0050670", "GO:0050671", "GO:0050870", "GO:0070661", "GO:0070663", "GO:0070665", "GO:1903039")

apd <- c("GO:0008064", "GO:0008154", "GO:0010639", "GO:0022411", "GO:0030042", "GO:0030832", "GO:0030833", "GO:0030834", "GO:0030835", "GO:0031334", "GO:0032271", "GO:0032507", "GO:0032984", "GO:0043242", "GO:0043244", "GO:0043254", "GO:0043624", "GO:0051261", "GO:0110053", "GO:1901879", "GO:1901880", "GO:1902903")

hsp <- c("GO:0002791", "GO:0002793", "GO:0030073", "GO:0032024", "GO:0035773", "GO:0035774", "GO:0045822", "GO:0046883", "GO:0046887", "GO:0050714", "GO:0050796", "GO:0051047", "GO:0051222", "GO:0061178", "GO:0071326", "GO:0071331", "GO:0090087", "GO:0090276", "GO:0090277", "GO:1903523", "GO:1904951")

wpc <- c("GO:0007596", "GO:0007599", "GO:0010543", "GO:0010544", "GO:0030193", "GO:0032102", "GO:0034109", "GO:0034110", "GO:0034111", "GO:0042060", "GO:0050817", "GO:0050818", "GO:0061041", "GO:0061045", "GO:0070527", "GO:0090330", "GO:0090331", "GO:1900046", "GO:1903034", "GO:1903035")
ntm <- c("GO:0001655", "GO:0001656", "GO:0001822", "GO:0001823", "GO:0030278", "GO:0032835", "GO:0060993", "GO:0072001", "GO:0072006", "GO:0072028", "GO:0072073", "GO:0072078", "GO:0072080", "GO:0072087", "GO:0072088", "GO:0072189", "GO:0090184")

scs <- c("GO:0000070", "GO:0000280", "GO:0000819", "GO:0007051", "GO:0007059", "GO:0007062", "GO:0007127", "GO:0045143", "GO:0048285", "GO:0051225", "GO:0051321", "GO:0061982", "GO:0071459", "GO:0098813", "GO:0140013", "GO:0140014", "GO:1902850")

tpi <- c("GO:0001508", "GO:0010644", "GO:0032410", "GO:0032413", "GO:0034763", "GO:0035637", "GO:0043266", "GO:0043267", "GO:0043271", "GO:0071805", "GO:0086001", "GO:0086002", "GO:0086011", "GO:0086064", "GO:1901379", "GO:1901380", "GO:1904063")

ppc <- c("GO:0009896", "GO:0010042", "GO:0010950", "GO:0010952", "GO:0031331", "GO:0032434", "GO:0042176", "GO:0045732", "GO:0045862", "GO:0061136", "GO:1900221", "GO:1903050", "GO:1903362", "GO:2000058")

mfc <- c("GO:0001539", "GO:0001578", "GO:0003341", "GO:0003351", "GO:0007271", "GO:0032222", "GO:0032224", "GO:0035082", "GO:0044782", "GO:0060271", "GO:0060285", "GO:0060294", "GO:0090660")

inn <- c("GO:0006913", "GO:0032386", "GO:0033157", "GO:0034504", "GO:0042306", "GO:0046822", "GO:0046824", "GO:0051169", "GO:1900180", "GO:1900182", "GO:1904589")
stc <- c("GO:0003014", "GO:0007200", "GO:0007204", "GO:0010459", "GO:0035025", "GO:0046579", "GO:0051057", "GO:0051281", "GO:0051480", "GO:0051482", "GO:0098801")

btr <- c("GO:0032091", "GO:0032092", "GO:0043392", "GO:0043393", "GO:0051098", "GO:0051100", "GO:0051101", "GO:2000677", "GO:2000678")
gpt <- c("GO:0000082", "GO:0044770", "GO:0044772", "GO:0044843", "GO:0045787", "GO:0051984", "GO:0090068", "GO:1901987", "GO:1902806")
rcc <- c("GO:0009166", "GO:0009261", "GO:0019439", "GO:0034655", "GO:0044270", "GO:0046434", "GO:0072523", "GO:1901292", "GO:1901361")

gim <- c("GO:0022612", "GO:0030879", "GO:0060056", "GO:0060443", "GO:0060444", "GO:0060603", "GO:0060644", "GO:0061180")
mts <- c("GO:0007519", "GO:0045844", "GO:0048636", "GO:0048643", "GO:0051149", "GO:0060538", "GO:1901861", "GO:1901863")
pch <- c("GO:0003298", "GO:0003301", "GO:0014743", "GO:0055021", "GO:0060420", "GO:0061049", "GO:0061050", "GO:0061051")

mdc <- c("GO:0007212", "GO:0071867", "GO:0071868", "GO:0071869", "GO:0071870", "GO:1903350", "GO:1903351")
nls <- c("GO:0009267", "GO:0009991", "GO:0031667", "GO:0031668", "GO:0031669", "GO:0042594", "GO:0071496")
rte <- c("GO:0003407", "GO:0046530", "GO:0048592", "GO:0048593", "GO:0060041", "GO:0060042", "GO:0090596")
vhd <- c("GO:0009615", "GO:0042509", "GO:0050691", "GO:0051607", "GO:0071359", "GO:0071360", "GO:0140546")

alm <- c("GO:0006301", "GO:0007611", "GO:0008306", "GO:0009314", "GO:0009416", "GO:0019985")
dos <- c("GO:0036474", "GO:1901031", "GO:1901032", "GO:1902883", "GO:1903202", "GO:1903206")
ipn <- c("GO:0050709", "GO:0051224", "GO:0090317", "GO:1903531", "GO:1903828", "GO:1904950")
nma <- c("GO:0014041","GO:0042551","GO:0045664","GO:0045665","GO:0048708","GO:1903429")

bto <- c("GO:0032200","GO:0033044","GO:0051052","GO:0071897","GO:2000278")
bpr <- c("GO:0031529","GO:0060491","GO:0097178","GO:0120032","GO:1900027")
csa <- c("GO:0090342","GO:0090343","GO:0090398","GO:2000772","GO:2000774")
igp <-c("GO:0001562","GO:0034341","GO:0042832","GO:0070673","GO:0071346")
mdh <- c("GO:0002573","GO:0030101","GO:0045637","GO:1902105","GO:1903706")
swc <- c("GO:0030431","GO:0042745","GO:0048520","GO:0050795","GO:0050802")

cws <- c("GO:0030111","GO:0030177","GO:0030178","GO:0060828")
cri <- c("GO:0002029","GO:0002031","GO:0008277","GO:0022401")
ecd <- c("GO:0042471","GO:0043583","GO:0048839","GO:0060113")
edt <- c("GO:0048557","GO:0048565","GO:0048566","GO:0055123")
epp <- c("GO:0019233","GO:0031644","GO:0098815","GO:2000463")
ees <- c("GO:0030198","GO:0043062","GO:0045229","GO:0071711")
mmc <- c("GO:0010753","GO:0019932","GO:0050848","GO:0050849")
sfr <- c("GO:0007338","GO:0009566","GO:2000241","GO:2000243")
smz <- c("GO:0001964","GO:0010043","GO:0071248","GO:0071294")
sam <- c("GO:0001952","GO:0007160","GO:0010810","GO:0048041")
sdh <- c("GO:0050803","GO:0050807","GO:0097284","GO:2000171")
sap <- c("GO:0001991","GO:0003081","GO:0044060","GO:0050886")
vsm <- c("GO:0048659","GO:0048660","GO:0048662","GO:1904706")

dsp <- c("GO:0007586","GO:0022600","GO:0050892")
emp <- c("GO:0007163","GO:0030010","GO:0090162")
jcm <- c("GO:0034331","GO:0043954","GO:0045216")
ntf <- c("GO:0001841","GO:0021915","GO:0042249")
pep <- c("GO:0006509","GO:0033619","GO:0042982")
spr <- c("GO:0060021","GO:0060022","GO:0062009")
scp <- c("GO:0006939","GO:0006940","GO:0014821")
sdg <- c("GO:0030307","GO:0045927","GO:0051155")
uan <- c("GO:0006914","GO:0044804","GO:0061919")

arm <- c("GO:0050854","GO:0050858")
ena <- c("GO:0048483","GO:0048484")
eec <- c("GO:0070371","GO:0070374")
esb <- c("GO:0030104","GO:0061436")
gpb <- c("GO:0010559","GO:1903018")
iml <- c("GO:0006690","GO:0043651")
moc <- c("GO:0007098","GO:0031023")
nrp <- c("GO:0006098","GO:0006740")
nsa <- c("GO:0010656","GO:0010664")
oep <- c("GO:0015748","GO:0015914")
pam <- c("GO:0007156","GO:0098742")
pmf <- c("GO:0000768","GO:0140253")
pga <- c("GO:0043087","GO:0043547")
plc <- c("GO:0007224","GO:0061512")
rts <- c("GO:0031929","GO:0032007")
rgc <- c("GO:0031960","GO:0051384")
spp <- c("GO:0006684","GO:0009395")
tmf <- c("GO:0021885","GO:0022029")
wtf <- c("GO:0006833","GO:0042044")




vertex_attr(amp)$AA <- case_when(vertex_attr(amp)$name %in% pim ~ "production immune mast",
                                 vertex_attr(amp)$name %in% brd ~ "break replication dna",
                                 vertex_attr(amp)$name %in% pcl ~ "proliferation cell selection",
                                 vertex_attr(amp)$name %in% apd ~ "actin polymerization depolymerization",
                                 vertex_attr(amp)$name %in% hsp ~ "hormone secretion peptide",
                                 vertex_attr(amp)$name %in% wpc ~ "wound platelet coagulation",
                                 vertex_attr(amp)$name %in% ntm ~ "nephron tubule morphogenesis",
                                 vertex_attr(amp)$name %in% scs ~ "sister chromatid segregation",
                                 vertex_attr(amp)$name %in% tpi ~ "transmembrane potassium ion",
                                 vertex_attr(amp)$name %in% ppc ~ "proteasomal proteolysis catabolic",
                                 vertex_attr(amp)$name %in% mfc ~ "movement fluid cilium",
                                 vertex_attr(amp)$name %in% inn ~ "intracellular nucleus nucleocytoplasmic",
                                 vertex_attr(amp)$name %in% stc ~ "signal transduction concentration",
                                 vertex_attr(amp)$name %in% btr ~ "binding transcription regulatory",
                                 vertex_attr(amp)$name %in% gpt ~ "g1 phase transition",
                                 vertex_attr(amp)$name %in% rcc ~ "ribonucleotide catabolic compound",
                                 vertex_attr(amp)$name %in% gim ~ "gland involution mammary",
                                 vertex_attr(amp)$name %in% mts ~ "muscle tissue skeletal",
                                 vertex_attr(amp)$name %in% pch ~ "physiological cardiac hypertrophy",
                                 vertex_attr(amp)$name %in% mdc ~ "monoamine dopamine catecholamine",
                                 vertex_attr(amp)$name %in% nls ~ "nutrient levels starvation",
                                 vertex_attr(amp)$name %in% rte ~ "retina type eye",
                                 vertex_attr(amp)$name %in% vhd ~ "virus host defense",
                                 vertex_attr(amp)$name %in% alm ~ "associative learning memory",
                                 vertex_attr(amp)$name %in% dos ~ "death oxygen species",
                                 vertex_attr(amp)$name %in% ipn ~ "intracellular protein negative",
                                 vertex_attr(amp)$name %in% nma ~ "neuron maturation astrocyte",
                                 vertex_attr(amp)$name %in% bto ~ "biosynthetic telomere organization",
                                 vertex_attr(amp)$name %in% bpr ~ "bounded projection ruffle",
                                 vertex_attr(amp)$name %in% csa ~ "cellular senescence aging",
                                 vertex_attr(amp)$name %in% igp ~ "interferon gamma protozoan",
                                 vertex_attr(amp)$name %in% mdh ~ "myeloid differentiation hemopoiesis",
                                 vertex_attr(amp)$name %in% swc ~ "sleep wake circadian",
                                 vertex_attr(amp)$name %in% cws ~ "canonical wnt signaling",
                                 vertex_attr(amp)$name %in% cri ~ "coupled receptor internalization",
                                 vertex_attr(amp)$name %in% ecd ~ "ear cell differentiation",
                                 vertex_attr(amp)$name %in% edt ~ "embryonic digestive tract",
                                 vertex_attr(amp)$name %in% epp ~ "excitatory postsynaptic potential",
                                 vertex_attr(amp)$name %in% ees ~ "external encapsulating structure",
                                 vertex_attr(amp)$name %in% mmc ~ "messenger mediated cgmp",
                                 vertex_attr(amp)$name %in% sfr ~ "single fertilization reproductive",
                                 vertex_attr(amp)$name %in% smz ~ "startle metal zinc",
                                 vertex_attr(amp)$name %in% sam ~ "substrate adhesion matrix",
                                 vertex_attr(amp)$name %in% sdh ~ "synapse dendrite hepatocyte",
                                 vertex_attr(amp)$name %in% sap ~ "systemic arterial pressure",
                                 vertex_attr(amp)$name %in% vsm ~ "vascular smooth muscle",
                                 vertex_attr(amp)$name %in% dsp ~ "digestive system process",
                                 vertex_attr(amp)$name %in% emp ~ "establishment maintenance polarity",
                                 vertex_attr(amp)$name %in% jcm ~ "junction component maintenance",
                                 vertex_attr(amp)$name %in% ntf ~ "neural tube formation",
                                 vertex_attr(amp)$name %in% pep ~ "precursor ectodomain proteolysis",
                                 vertex_attr(amp)$name %in% spr ~ "secondary palate roof",
                                 vertex_attr(amp)$name %in% scp ~ "smooth contraction phasic",
                                 vertex_attr(amp)$name %in% sdg ~ "striated differentiation growth",
                                 vertex_attr(amp)$name %in% uan ~ "utilizing autophagy nucleus",
                                 vertex_attr(amp)$name %in% arm ~ "antigen receptor mediated",
                                 vertex_attr(amp)$name %in% ena ~ "enteric nervous autonomic",
                                 vertex_attr(amp)$name %in% eec ~ "erk1 erk2 cascade",
                                 vertex_attr(amp)$name %in% esb ~ "establishment skin barrier",
                                 vertex_attr(amp)$name %in% gpb ~ "glycoprotein process biosynthetic",
                                 vertex_attr(amp)$name %in%  iml ~ "icosanoid metabolic linoleic",
                                 vertex_attr(amp)$name %in% moc ~ "microtubule organizing center",
                                 vertex_attr(amp)$name %in% nrp ~ "nadph regeneration pentose",
                                 vertex_attr(amp)$name %in% nsa ~ "negative striated apoptotic",
                                 vertex_attr(amp)$name %in% oep ~ "organophosphate ester phospholipid",
                                 vertex_attr(amp)$name %in% pam ~ "plasma adhesion molecules",
                                 vertex_attr(amp)$name %in% pmf ~ "plasma membrane fusion",
                                 vertex_attr(amp)$name %in% pga ~ "positive gtpase activity",
                                 vertex_attr(amp)$name %in% plc ~ "protein localization cilium",
                                 vertex_attr(amp)$name %in% rts ~ "regulation tor signaling",
                                 vertex_attr(amp)$name %in% rgc ~ "response glucocorticoid corticosteroid",
                                 vertex_attr(amp)$name %in% spp ~ "sphingomyelin process phospholipid",
                                 vertex_attr(amp)$name %in% tmf ~ "telencephalon migration forebrain",
                                 vertex_attr(amp)$name %in% wtf ~ "water transport fluid")

long_list_to_fill_amp <- lapply(seq_along(unique(vertex_attr(amp)$AA)),function(i) i)

AA_Clusters <- lapply(seq_along(vertex_attr(amp)$AA),function(n) {
  #unj_2 <- list()
  #print(n)
  if_else(is.na(vertex_attr(amp)$AA[[n]]) == TRUE,
          long_list_to_fill_amp[n] <- vertex_attr(amp)$`EnrichmentMap::GS_DESCR`[[n]],
          long_list_to_fill_amp[n] <- vertex_attr(amp)$AA[[n]])
})


vertex_attr(amp)$Clusters <- unlist(AA_Clusters)


#################################
#   Uninjured                   #
#################################
unj <- read_graph(file = "./data/cleaned_unique_go_bp_v2__Uninjured_vs_Sham.graphml",
                  format = "graphml")


pdp <- c("GO:0046939","GO:0046365","GO:0019320","GO:0009185","GO:0009179","GO:0009135","GO:0009132","GO:0006757","GO:0006096")
btd <- c("GO:0110148", "GO:0051216", "GO:0048286", "GO:0045778", "GO:0031214", "GO:0030282", "GO:0002062")
eap <- c("GO:2001237", "GO:2001236", "GO:0097191", "GO:0070301")
afa <- c("GO:0030038","GO:0051492","GO:0043149")
egr <- c("GO:1901185","GO:0042059","GO:0007175")
fss <- c("GO:0071498","GO:0034405")
nkc <- c("GO:0007566","GO:0007565")
vertex_attr(unj)$AA <- vertex_attr(unj)$`EnrichmentMap::GS_DESCR`
vertex_attr(unj)$AA <- case_when(vertex_attr(unj)$name %in% pdp ~ "purine diphosphate process",
                                 vertex_attr(unj)$name %in% btd ~ "biomineral tissue development",
                                 vertex_attr(unj)$name %in% eap ~ "extrinsic apoptotic pathway",
                                 vertex_attr(unj)$name %in% afa ~ "actin filament assembly",
                                 vertex_attr(unj)$name %in% egr ~ "epidermal growth negative",
                                 vertex_attr(unj)$name %in% fss ~ "fluid shear stress",
                                 vertex_attr(unj)$name %in% nkc ~ "natural killer cell")

long_list_to_fill_unj <- lapply(seq_along(unique(vertex_attr(unj)$AA)),function(i) i)
AA_Clusters <- lapply(seq_along(vertex_attr(unj)$AA),function(n) {
  if_else(is.na(vertex_attr(unj)$AA[[n]]) == TRUE,
          long_list_to_fill_unj[n] <- vertex_attr(unj)$`EnrichmentMap::GS_DESCR`[[n]],
          long_list_to_fill_unj[n] <- vertex_attr(unj)$AA[[n]])
})

vertex_attr(unj)$Clusters <- unlist(AA_Clusters)



#################################
#     Core regeneration         #
#################################
core <- read_graph(file = "./data/heart_regeneration_core_mm_fil_cleaned_Abl_Sham__gobp.graphml",
                   format = "graphml")

# vertex_attr(core)$name
# vertex_attr(core)$`EnrichmentMap::GS_DESCR`
cmc <- c("GO:0000075", "GO:0007088", "GO:0007091", "GO:0007093", "GO:0007094", "GO:0007096", "GO:0010639", "GO:0010948", "GO:0010965", "GO:0030071", "GO:0031577", "GO:0033044", "GO:0033045", "GO:0033046", "GO:0033047", "GO:0033048", "GO:0044770", "GO:0044772", "GO:0044784", "GO:0045786", "GO:0045839", "GO:0045841", "GO:0051304", "GO:0051306", "GO:0051783", "GO:0051784", "GO:0051983", "GO:0051985", "GO:0071173", "GO:0071174", "GO:0090231", "GO:0090266", "GO:1901976", "GO:1901987", "GO:1901988", "GO:1901990", "GO:1901991", "GO:1902099", "GO:1902100", "GO:1903504", "GO:1905818", "GO:1905819", "GO:2000816", "GO:2001251")
msm <- c("GO:1903046", "GO:0140014", "GO:0140013", "GO:0098813", "GO:0090306", "GO:0070601", "GO:0070192", "GO:0061983", "GO:0051988", "GO:0051321", "GO:0051310", "GO:0051303", "GO:0051177", "GO:0050000", "GO:0048285", "GO:0045144", "GO:0031145", "GO:0008608", "GO:0007135", "GO:0007062", "GO:0007059", "GO:0000819", "GO:0000280", "GO:0000212", "GO:0000070")
wcp <- c("GO:0198738", "GO:0090263", "GO:0090090", "GO:0060828", "GO:0060070", "GO:0035567", "GO:0030178", "GO:0030177", "GO:0030111", "GO:0016055")
cpp <- c("GO:2000060", "GO:1903364", "GO:1903052", "GO:1901800", "GO:0045732", "GO:0042176", "GO:0032436", "GO:0032434", "GO:0009896")
bom <- c("GO:0110148", "GO:0031214", "GO:0030282", "GO:0030278", "GO:0001649", "GO:0001503")
bbc <- c("GO:1903036","GO:0050878","GO:0050817","GO:0042060","GO:0007599","GO:0007596")
dte <- c("GO:0150063","GO:0090596","GO:0060041","GO:0048880","GO:0043010","GO:0001654")
mas <- c("GO:0045785","GO:0031589","GO:0010812","GO:0010811","GO:0010810","GO:0001953")
aab <- c("GO:1901607","GO:0046394","GO:0016053","GO:0009070","GO:0008652")
ccc <- c("GO:0061448","GO:0061036","GO:0061035","GO:0051216","GO:0002062")
csg <- c("GO:0048672","GO:0048639","GO:0045927","GO:0030307","GO:0010720")
ebs <- c("GO:1904888","GO:0060350","GO:0060349","GO:0060348","GO:0048705")
sao <- c("GO:1902850","GO:0090307","GO:0051225","GO:0007052","GO:0007051")
tma <- c("GO:0090132","GO:0090130","GO:0010763","GO:0010631","GO:0001667")
bvr <- c("GO:0048771","GO:0046849","GO:0034103","GO:0001974")
emd <- c("GO:0060485","GO:0048762","GO:0010718","GO:0010717")
ese <- c("GO:0045229","GO:0043062","GO:0030199","GO:0030198")
gft <- c("GO:1903844","GO:0090287","GO:0040037","GO:0017015")
mcc <- c("GO:0061640","GO:0051302","GO:0000910","GO:0000281")
pcc <- c("GO:0050766","GO:0050764","GO:0043277","GO:0006909")
pdr <- c("GO:0048008","GO:0035791","GO:0010642","GO:0010640")
smm <- c("GO:0014911","GO:0014910","GO:0014909","GO:0014812")
adi <- c("GO:2001233","GO:0097193","GO:0008630")
ama <- c("GO:0060840","GO:0048844","GO:0035904")
crg <- c("GO:0071622","GO:0071621","GO:0060326")
die <- c("GO:1990138","GO:0060560","GO:0048588")
fgm <- c("GO:0048332","GO:0001707","GO:0001704")
gaa <- c("GO:0097485","GO:0007411","GO:0007409")
lhd <- c("GO:0097421","GO:0061008","GO:0001889")
mdm <- c("GO:1904181","GO:0051901","GO:0051900")
smp <- c("GO:0048660","GO:0048659","GO:0033002")
ser <- c("GO:0035815","GO:0007588","GO:0003014")
tka <- c("GO:0071900","GO:0045860","GO:0033674")
uat <- c("GO:1904666","GO:0051443","GO:0051438")
afb <- c("GO:0032970","GO:0032956")
abc <- c("GO:1900221","GO:0097242")
cjd <- c("GO:0045216","GO:0002934")
cccc <- c("GO:0010457","GO:0007098")
cccp <- c("GO:2000343","GO:0001819")
fgg <- c("GO:0048477","GO:0007292")
gpc <- c("GO:0044839","GO:0000086")
mcs <- c("GO:0051299","GO:0007100")
nlb <- c("GO:0051055","GO:0046890")
nppp <- c("GO:2000179","GO:2000177")
npr <- c("GO:0031102","GO:0031099")
ocl <- c("GO:0061842","GO:0051642")
psp <- c("GO:0018209","GO:0018105")
pcp <- c("GO:1901989","GO:0090068")
ros <- c("GO:0072593","GO:0006801")
rpb <- c("GO:0051098","GO:0043393")
sbp <- c("GO:0008217","GO:0003073")






vertex_attr(core)$AA <- vertex_attr(core)$`EnrichmentMap::GS_DESCR`
vertex_attr(core)$AA <- case_when(vertex_attr(core)$name %in% cmc ~ "checkpoint mitotic cycle",
                                  vertex_attr(core)$name %in% msm ~ "meiosis sister meiotic",
                                  vertex_attr(core)$name %in% wcp ~ "wnt canonical pathway",
                                  vertex_attr(core)$name %in% cpp ~ "catabolic proteasomal protein",
                                  vertex_attr(core)$name %in% bom ~ "biomineral ossification mineralization",
                                  vertex_attr(core)$name %in% bbc ~ "body blood coagulation",
                                  vertex_attr(core)$name %in% dte ~ "development type eye",
                                  vertex_attr(core)$name %in% mas ~ "matrix adhesion substrate",
                                  vertex_attr(core)$name %in% aab ~ "acid amino biosynthetic",
                                  vertex_attr(core)$name %in% ccc ~ "cartilage connective chondrocyte",
                                  vertex_attr(core)$name %in% csg ~ "collateral sprouting growth",
                                  vertex_attr(core)$name %in% ebs ~ "endochondral bone skeletal",
                                  vertex_attr(core)$name %in% sao ~ "spindle assembly organization",
                                  vertex_attr(core)$name %in% tma ~ "tissue migration ameboidal",
                                  vertex_attr(core)$name %in% bvr ~ "blood vessel remodeling",
                                  vertex_attr(core)$name %in% emd ~ "epithelial mesenchymal differentiation",
                                  vertex_attr(core)$name %in% ese ~ "external structure extracellular",
                                  vertex_attr(core)$name %in% gft ~ "growth factor transforming",
                                  vertex_attr(core)$name %in% mcc ~ "mitotic cytokinesis cytoskeleton",
                                  vertex_attr(core)$name %in% pcc ~ "phagocytosis cell clearance",
                                  vertex_attr(core)$name %in% pdr ~ "platelet derived receptor",
                                  vertex_attr(core)$name %in% smm ~ "smooth muscle migration",
                                  vertex_attr(core)$name %in% adi ~ "apoptotic damage intrinsic",
                                  vertex_attr(core)$name %in% ama ~ "artery morphogenesis aorta",
                                  vertex_attr(core)$name %in% crg ~ "chemotaxis regulation granulocyte",
                                  vertex_attr(core)$name %in% die~ "developmental involved extension",
                                  vertex_attr(core)$name %in% fgm~ "formation germ mesoderm",
                                  vertex_attr(core)$name %in% gaa ~ "guidance axon axonogenesis",
                                  vertex_attr(core)$name %in% lhd ~ "liver hepaticobiliary development",
                                  vertex_attr(core)$name %in% mdm ~ "membrane depolarization mitochondrial",
                                  vertex_attr(core)$name %in% smp ~ "smooth muscle proliferation",
                                  vertex_attr(core)$name %in% ser ~ "sodium excretion renal",
                                  vertex_attr(core)$name %in% tka ~ "threonine kinase activity",
                                  vertex_attr(core)$name %in% uat ~ "ubiquitin activity transferase",
                                  vertex_attr(core)$name %in% afb ~ "actin filament based",
                                  vertex_attr(core)$name %in% abc ~ "amyloid beta clearance",
                                  vertex_attr(core)$name %in% cjd  ~ "cell junction desmosome",
                                  vertex_attr(core)$name %in% cccc ~ "centriole cohesion centrosome",
                                  vertex_attr(core)$name %in% cccp ~ "chemokine cytokine production",
                                  vertex_attr(core)$name %in% fgg  ~ "female gamete generation",
                                  vertex_attr(core)$name %in% gpc ~ "g2 phase cell",
                                  vertex_attr(core)$name %in% mcs ~ "mitotic centrosome separation",
                                  vertex_attr(core)$name %in% nlb ~ "negative lipid biosynthetic",
                                  vertex_attr(core)$name %in% nppp ~ "neural precursor proliferation",
                                  vertex_attr(core)$name %in% npr ~ "neuron projection regeneration",
                                  vertex_attr(core)$name %in% ocl ~ "organizing center localization",
                                  vertex_attr(core)$name %in% psp ~ "peptidyl serine phosphorylation",
                                  vertex_attr(core)$name %in% pcp ~ "positive cycle phase",
                                  vertex_attr(core)$name %in% ros ~ "reactive oxygen species",
                                  vertex_attr(core)$name %in% rpb ~ "regulation protein binding",
                                  vertex_attr(core)$name %in% sbp ~ "systemic blood pressure")





long_list_to_fill_core <- lapply(seq_along(unique(vertex_attr(core)$AA)),function(i) i)

AA_Clusters <- lapply(seq_along(vertex_attr(core)$AA),function(n) {
  if_else(is.na(vertex_attr(core)$AA[[n]]) == TRUE,
          long_list_to_fill_core[n] <- vertex_attr(core)$`EnrichmentMap::GS_DESCR`[[n]],
          long_list_to_fill_core[n] <- vertex_attr(core)$AA[[n]])
})
vertex_attr(core)$Clusters <- unlist(AA_Clusters)


###############################
#       Load counts           #
###############################

unj_counts <- read.delim(file = "./data/cleaned_unique_go_bp_v2__Uninjured_vs_Sham_normalized_counts_expression_symbol_and_condition_enmap.txt",header = TRUE,sep = "\t")
unj_counts[nrow(unj_counts)+1,] <- 0
unj_counts[nrow(unj_counts)+1,] <- 0
unj_counts[nrow(unj_counts),"Symbols"] <- "Select_a_node"
unj_counts[nrow(unj_counts)-1,"Symbols"] <- "Select_a_node"


abl_counts <- read.delim(file = "./data/cleaned_unique_go_bp_v2__Ablated_vs_Sham_normalized_counts_expression_symbol_and_condition_enmap.txt",header = TRUE,sep = "\t")
abl_counts[nrow(abl_counts)+1,] <- 0
abl_counts[nrow(abl_counts)+1,] <- 0
abl_counts[nrow(abl_counts),"Symbols"] <- "Select_a_node"
abl_counts[nrow(abl_counts)-1,"Symbols"] <- "Select_a_node"

amp_counts <- read.delim(file = "./data/cleaned_unique_go_bp_v2__Resection_vs_Sham_normalized_counts_expression_symbol_and_condition_enmap.txt",header = TRUE,sep = "\t")
amp_counts[nrow(amp_counts)+1,] <- 0
amp_counts[nrow(amp_counts)+1,] <- 0
amp_counts[nrow(amp_counts),"Symbols"] <- "Select_a_node"
amp_counts[nrow(amp_counts)-1,"Symbols"] <- "Select_a_node"

cry_counts <- read.delim(file = "./data/Cryonjury_vs_Sham_mouse_enriched_all_clean.txt",header = TRUE,sep = "\t")
cry_counts[nrow(cry_counts)+1,] <- 0
cry_counts[nrow(cry_counts)+1,] <- 0
cry_counts[nrow(cry_counts),"Symbols"] <- "Select_a_node"
cry_counts[nrow(cry_counts)-1,"Symbols"] <- "Select_a_node"


s4c <- read.delim(file = "./data/colData_s4c.txt",header = TRUE,sep = "\t")


#sham_counts <- cry_counts[,c(rownames(s4c %>% filter(Condition=="Sham")),"MGI_Symbol")]



core_counts <- read.delim("./data/core_selected_counts_148_Table.txt",sep = "\t",header = TRUE)
core_counts[nrow(core_counts)+1,] <- 0
core_counts[nrow(core_counts)+1,] <- 0
core_counts[nrow(core_counts),"MGI_Symbol"] <- "Select_a_node"
core_counts[nrow(core_counts)-1,"MGI_Symbol"] <- "Select_a_node"


enrich_go_bo_df <- read.delim("./data/Mmus_clean_gobp_enrichment_v2_abl_sham.txt",sep = "\t",header = TRUE)
enrich_go_bo_df$Group[enrich_go_bo_df$Group=="Amputation" ] <- "Resection"


core_reg <- read.delim(file = "./data/heart_regeneration_core_mm_fil_cleaned_Abl_Sham__gobp.txt",sep = "\t",header = TRUE)


#gene seeker
abl_fil <- read.delim(file = "./data/Ablated_vs_Sham_mouse_enriched_fil_clean.txt",header = TRUE,sep = "\t")
amp_fil <- read.delim(file = "./data/Resection_vs_Sham_mouse_enriched_fil_clean.txt",header = TRUE,sep = "\t")
cryo_fil <- read.delim(file = "./data/Cryonjury_vs_Sham_mouse_enriched_fil_clean.txt",header = TRUE,sep = "\t")
unj_fil <- read.delim(file = "./data/Uninjured_vs_Sham_mouse_enriched_fil_clean.txt",header = TRUE,sep = "\t")

unj_fil_f <- unj_fil |> dplyr::filter(!is.na(MGI_Symbol)) |> dplyr::select(MGI_Symbol,starts_with("SRR"),Description)
amp_fil_f <- amp_fil |> dplyr::filter(!is.na(MGI_Symbol)) |> dplyr::select(MGI_Symbol,starts_with("SRR"),Description)
abl_fil_f <- abl_fil |> dplyr::filter(!is.na(MGI_Symbol)) |> dplyr::select(MGI_Symbol,starts_with("SRR"),Description)
cryo_fil_f <- cryo_fil |> dplyr::filter(!is.na(MGI_Symbol)) |> dplyr::select(MGI_Symbol,starts_with("SRR"),Description)

gene_seeker_table <- data.frame()
gene_seeker_table <- rbind(gene_seeker_table,abl_fil_f,amp_fil_f,unj_fil_f,cryo_fil_f)

#choices = c(as.character(sort(unique(gene_seeker_table$MGI_Symbol)))),


gene_seeker_table <- gene_seeker_table |> dplyr::group_by(MGI_Symbol) |> distinct(MGI_Symbol,.keep_all = TRUE)




#log2fc heatmap of the core regeneration genes
files <- list.files(all.files = TRUE,
                    path = "./data",
                    #path = "D:/PhD/Projects/prsa/Outputs/r/Cytoscape/Cleaned_Annotation_Ensembl_v4/dr2mmus/",
                    pattern = "fil_clean.txt",
                    recursive = TRUE, 
                    full.names = TRUE)

tables_counts_stats <- lapply(files,function(n) read.delim(n,header = TRUE,sep = "\t"))

name_files <- gsub(x = files,
                   pattern = c("./data/"),
                   replacement = "")
name_files <- gsub(x = name_files,
                   pattern = c("_enriched_fil_clean.txt"),
                   replacement = "")

names(tables_counts_stats) <- name_files

#create df
cryo_df <- tables_counts_stats$Cryonjury_vs_Sham_mouse %>% dplyr::filter(MGI_Symbol %in% core_counts$MGI_Symbol)
cryo_df <- cryo_df %>% dplyr::select(MGI_Symbol,log2FoldChange,starts_with("SRR"))
#dim(cryo_df)


amp_df <- tables_counts_stats$Resection_vs_Sham_mouse %>% filter(MGI_Symbol %in% core_counts$MGI_Symbol)
amp_df <- amp_df %>% dplyr::select(MGI_Symbol,log2FoldChange,starts_with("SRR"))
# print(dim(amp_df))

abl_df <- tables_counts_stats$Ablated_vs_Sham_mouse %>% filter(MGI_Symbol %in% core_counts$MGI_Symbol)
abl_df <- abl_df %>% dplyr::select(MGI_Symbol,log2FoldChange,starts_with("SRR"))
#dim(abl_df)


#abl_df %>% group_by(MGI_Symbol) %>% count(MGI_Symbol)  %>% filter(n > 1)
#abl_df %>% group_by(MGI_Symbol) %>% distinct(MGI_Symbol) %>% count(MGI_Symbol) %>% filter(n > 1)

abl_df_2 <- abl_df %>% group_by(MGI_Symbol) %>% distinct(MGI_Symbol,.keep_all = TRUE) %>% dplyr::arrange(MGI_Symbol)
amp_df_2 <- amp_df %>% group_by(MGI_Symbol) %>% distinct(MGI_Symbol,.keep_all = TRUE) %>% dplyr::arrange(MGI_Symbol)
cryo_df_2 <- cryo_df %>% group_by(MGI_Symbol) %>% distinct(MGI_Symbol,.keep_all = TRUE) %>% dplyr::arrange(MGI_Symbol)


#rm(counts_148_selected_without_doubles)
all.equal(cryo_df_2$MGI_Symbol,abl_df_2$MGI_Symbol)
counts_148_selected_without_doubles <- as.data.frame(amp_df_2$log2FoldChange)
counts_148_selected_without_doubles <- cbind.data.frame(counts_148_selected_without_doubles,abl_df_2$log2FoldChange)
counts_148_selected_without_doubles <- cbind.data.frame(counts_148_selected_without_doubles,cryo_df_2$log2FoldChange)

rownames(counts_148_selected_without_doubles) <- amp_df_2$MGI_Symbol
colnames(counts_148_selected_without_doubles) <- c("Resection","Ablation","Cryoinjury")


ui <- secure_app(head_auth = tags$script(inactivity),
                 dashboardPage(skin = "purple",
                               dashboardHeader(title = "Heart Injuries"),
                               dashboardSidebar(width = 300,
                                                sidebarMenu(
                                                  menuItem("Home",tabName = "home",icon = icon("home")),
                                                  menuItem("Resection vs Sham",tabName = "resected", icon = icon("scissors",lib = "glyphicon")),#,icon = icon("calendar",lib = "glyphicon")),
                                                  menuItem("Ablation vs Sham",tabName = "ablated",icon = icon("tag",lib = "glyphicon")),
                                                  menuItem("Uninjured vs Sham",tabName = "uninjured",icon = icon("ok",lib = "glyphicon")),
                                                  menuItem("Heart Core Regeneration Common Genes ",tabName = "heartcore",icon = icon("heart-empty",lib = "glyphicon")),
                                                  menuItem("PubMed Query",tabName = "kuk",icon = icon("list",lib = "glyphicon")),
                                                  menuItem("Core Regeneration Log2FC",tabName = "core_l2fc",icon = icon("bookmark",lib = "glyphicon")),
                                                  menuItem("Genes Seeker",tabName = "gs",icon = icon("search",lib = "glyphicon")),
                                                  menuItem("About",tabName = "about",icon = icon("thumbtack"))
                                                )),
                               dashboardBody(
                                 tabItems(
                                   tabItem(tabName = "home",
                                           #home section and markdown
                                           includeMarkdown("./www/markdown_home.md")),
                                   # #next tab
                                   
                                   #https://stackoverflow.com/questions/70689513/how-to-have-shiny-dashboard-box-fill-the-entire-width-of-my-dashbaord
                                   #next tab
                                   tabItem(tabName = "resected",
                                           fluidRow(
                                             column(width = 12, textOutput("verb_amp")),#tags$head(tags$style("#text1{color: red;font-size: 20px;font-style: italic;}")) %>% withSpinner()),
                                             column(width = 12, visNetworkOutput(outputId = "resected_net",width = "auto", height = "600px") %>% withSpinner()),
                                             column(width = 12, dataTableOutput(outputId = "goTable_amp") %>% withSpinner()),
                                             column(width = 12,valueBoxOutput(width = "100%",outputId = "genes_go_amp") %>% withSpinner()),
                                             column(width = 12,plotlyOutput(outputId = "amphm",width = "100%",height = "100%") %>% withSpinner())
                                             
                                           )
                                   ),
                                   #next tab
                                   tabItem(tabName = "ablated",
                                           fluidRow(
                                             column(width = 12, textOutput("verb_abl")),#tags$head(tags$style("#text1{color: red;font-size: 20px;font-style: italic;}")) %>% withSpinner()),
                                             column(width = 12, visNetworkOutput(outputId = "ablated_net",width = "auto",height = "600px") %>% withSpinner()),
                                             column(width = 12, dataTableOutput(outputId = "goTable_abl") %>% withSpinner()),
                                             column(width = 12,valueBoxOutput(width = 12,outputId = "genes_go_abl") %>% withSpinner()),
                                             column(width = 12,plotlyOutput(outputId = "ablhm",width = "100%",height = "100%") %>% withSpinner())
                                             
                                           )
                                   ),
                                   #next tab
                                   tabItem(tabName = "uninjured",
                                           fluidRow(
                                             column(width = 12, textOutput("verb_unj")),#tags$head(tags$style("#text1{color: red;font-size: 20px;font-style: italic;}")) %>% withSpinner()),
                                             column(width = 12, visNetworkOutput(outputId = "uninjured_net",width = "auto",height = "600px") %>% withSpinner()),
                                             column(width = 12, dataTableOutput(outputId = "goTable_unj") %>% withSpinner()),
                                             column(width = 12,valueBoxOutput(width = 12,outputId = "genes_go_unj") %>% withSpinner()),
                                             column(width = 12,plotlyOutput(outputId = "unjhm",width = "100%",height = "100%") %>% withSpinner())
                                             
                                           )
                                   ),
                                   #next tab
                                   tabItem(tabName = "heartcore",
                                           fluidRow(
                                             column(width = 12, textOutput("verb_core") %>% withSpinner()),
                                             column(width = 12, visNetworkOutput(outputId = "heartcore_net",width = "auto",height = "600px") %>% withSpinner()),
                                             column(width = 12, dataTableOutput(outputId = "goTable_core") %>% withSpinner()),
                                             column(width = 12,valueBoxOutput(width = "auto",outputId = "genes_go_core") %>% withSpinner()),
                                             column(width = 12,plotlyOutput(outputId = "hchm",width = "100%",height = "100%") %>% withSpinner())
                                             
                                           )
                                   ),
                                   #next tab
                                   tabItem(tabName = "kuk",
                                           fluidRow(
                                             column(width = 12,plotlyOutput(outputId = "known", width = "100%", height = "1024px") %>% withSpinner()),
                                             column(width = 12,plotlyOutput(outputId = "unknown", width = "100%", height = "1024px") %>% withSpinner())         
                                           )
                                   ),
                                   #next tab
                                   tabItem(tabName = "core_l2fc",
                                           fluidRow(
                                             column(width = 12,plotlyOutput(outputId = "core_l2fc_hm", width = "100%", height = "2048px") %>% withSpinner())         
                                           )
                                   ),
                                   #next tab
                                   tabItem(tabName = "gs",
                                           fluidRow(
                                             column(width=4,selectizeInput(inputId = "gene",
                                                                           choices = NULL,
                                                                           label = "Select at least 2 Gene Symbol",
                                                                           #choices = c(as.character(sort(unique(gene_seeker_table$MGI_Symbol)))),
                                                                           selected = c("Foxm1"),
                                                                           multiple = TRUE)),
                                             column(width = 12, dataTableOutput(outputId = "geneseeker_table") %>% withSpinner()),
                                             column(width = 12,plotlyOutput(outputId = "geneseeker_plot",width = "100%",height = "100%") %>% withSpinner()),
                                           )
                                   ),
                                   #next tab
                                   tabItem(tabName = "about",
                                           h2("Created and scratched by: "),
                                           
                                           fluidRow(
                                             
                                             fluidRow(
                                               column(width = 8,
                                                      h3(tags$a("Marius Botos",
                                                                href = "https://sporgelum.github.io/")),
                                                      imageOutput("marius")
                                               ),
                                               column(width = 4,imageOutput("img7")),
                                               style='padding-left:20px; padding-right:0px; padding-top:0px; padding-bottom:25px')
                                           ),
                                           
                                           fluidRow(
                                             column(width = 8,
                                                    h3(tags$a("Prateek Arora",
                                                              href = "https://www.anatomie.unibe.ch/ueber_uns/team/detail/index_ger.php?id=578")),
                                                    imageOutput("prateek"),
                                                    style='padding-left:20px; padding-right:0px; padding-top:0px; padding-bottom:25px'),
                                             column(width = 4,imageOutput("pa1"))),
                                           
                                           fluidRow(
                                             column(width = 8,
                                                    h3(tags$a("Panagiotis Chouvardas",
                                                              href = "https://www.urogenus-research.org/team/chouvardas-panagiotis")),
                                                    imageOutput("panos"),
                                                    style='padding-left:20px; padding-right:0px; padding-top:0px; padding-bottom:25px'),
                                             column(width = 2,imageOutput("pc1")),
                                             column(width = 2,imageOutput("pc2"))),
                                           
                                           fluidRow(
                                             column(width = 8,
                                                    h3(tags$a("Nadia Mercarder",
                                                              href = "https://www.anatomie.unibe.ch/about_us/management/detail/index_eng.php?id=449")),
                                                    imageOutput("nadia"),
                                                    style='padding-left:20px; padding-right:0px; padding-top:0px; padding-bottom:25px'),
                                             column(width = 2,imageOutput("nd1")),
                                             column(width = 2,imageOutput("nd2"))
                                             
                                           )
                                   )
                                 )# end tabItems
                                 # Add download button 
                               )# end dashboardBody
                 )# end dashboardPage
)#end all
