# Seurat - Guided Clustering Tutorial

library(dplyr)
library(Seurat)
library(patchwork)

# ------------------------------------------------------------------------------
#  PBMC 데이터셋 불러오기 
pbmc.data <- Read10X(data.dir = "~/Desktop/filtered_gene_bc_matrices/hg19/")

# ------------------------------------------------------------------------------
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

# ------------------------------------------------------------------------------
# pbmc 살펴보기 
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

#CD3D  4 . 10 . . 1 2 3 1 . . 2 7 1 . . 1 3 . 2  3 . . . . . 3 4 1 5
#TCL1A . .  . . . . . . 1 . . . . . . . . . . .  . 1 . . . . . . . .
#MS4A1 . 6  . . . . . . 1 1 1 . . . . . . . . . 36 1 2 . . 2 . . . .

# . : 분자가 발견되지 않았음 , scRNA-seq 대부분은 0이다. 

# ------------------------------------------------------------------------------
dense.size <- object.size(as.matrix(pbmc.data))
dense.size

sparse.size <- object.size(pbmc.data)
sparse.size

# dense.size vs sparse.size
# -> 같은 데이터를 다르게 표현했을 때 (희소, 밀집) 메모리를 얼마나 차지하는지 비교하는 코드
# -> dense.size : dense matrix -> 메모리 많이 사용 
# -> sparse.size : sparse matrix -> 메모리 절약 

# ------------------------------------------------------------------------------
### Standard pre-processing workflow : single cell RNA 시퀀싱 데이터를 전처리하는 과정 

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc
# PercentageFeatureSet : 지정한 유전자들의 발현량이 전체에서 차지하는 비율 계산
#  pattern = "^MT-" -> MT-로 시작하는 유전자들만 고른다. 
# pbmc[["percent.mt]] -> 계산된 결과를 percent.mt라는 이름으로 메타데이터에 저장해줌 

# -> 각 세포마다 미토콘드리아 유전자 발현 비율을 계산해서 붙여주는 것 ! 

# print 해보쟈  / Q 미코 콘드리아 비율이 많다고 판단하는 기준?? 
pbmc[["percent.mt"]]


head(pbmc@meta.data, 5)

# 세포 필터링 단계

# 시각화를 해보쟈
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Feature간의 관계를 알아볼 때 FeatureScatter 사용 !  
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt") 
plot2 <- FeatureScatter(pbmc, feature1= "nCount_RNA", feature2 = "nFeature_RNA")

plot1 + plot2

pbmc <- subset(pbmc,subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

### 정규화 
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# pbmc <- NormalizeData(pbmc) 이렇게 해도 똑같음 
# 이 방식은 모든 세포가 가진 RNA 총량이 같다고 가정한다. 

# SCTrasform() 이것도 사용해보기 

# ----
###  Identification of highly variable features (feature selection)
# 변화가 많은 유전자를 고르는 단계
# 세포마다 발현량 차이가 큰 유전자들을 뽑아내는 과정 

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst",nfeatures = 2000)

# Identify the 10 most highly variable genes

top10 <- head(VariableFeatures(pbmc),10)

# plot variable features with and without labels 

plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


# ----
###  Scaling the data : 차원 축소 분석(PCA) 를 하기 전에 필요한 단계

all.genes  <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes) # 전체 유전자에 대해 평균0, 분산1이 되도록 규모를 맞춰줌 

# 또한 원치 않은 변동 요인(unwanted variation) 를 제거(regress out)할 때도 사용할 수 있다. 
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")

# ----
### Perform linear dimensional reduction
# scRNA-seq 데이터에서 PCA(주성분분석)를 수행하고 결과를 시각화하는 과정

# 1. PCA란 ? 
# - 여러 유전자의 발현값(수천 ~ 만개)을 몇개의 주성분(PC)라는 축으로 줄여서 표현하는 방법
# - 이 주성분들은 데이터 안에 있는 가장 큰 변동성을 담고 있다.
# - 복잡한 데이터를 단순하게 요약해주는 차원 축소 기법

# pca 수행, 
# features 인자로는 분석할 유전자를 지정함, 보통 변동성 큰 유전자들 
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# pca 결과 출력
print(pbmc[["pca"]], dims= 1:5, nfeatures=5)

# 시각화 
# 그래프 보는 법 ?? 

# PC1과 PC2에 가장 큰 영향을 주는 유전자들을 막대 그래프로 보여준다. 
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca") 

# 세포들을 PC1과 PC2축에 점을 찍어 보여주는 산점도 
DimPlot(pbmc,reduction = "pca") + NoLegend()

# PC1에 따른 세포와 유전자들의 발현 태펀을 히트맵으로 보여줌 
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# ----
### Determine the ‘dimensionality’ of the dataset 
# PCA 이후에 몇개의 주성분을 사용할 지 결정하는 과정 
ElbowPlot(pbmc)

# - 처음에는 ElbowPlot에서 꺾이는 지점 + 조금 더 늘려서 쓰는게 좋다.
# - 10 ~ 30개 정도면 대부분의 분석에서 잘 작동한다. 

# ----
### Cluster the cells 
# Seurat에서 세포들을 클러스터링(군집화)하는 과정 
# - graph based clustering 을 사용한다. 
# - 먼저 유사한 세포들끼지 k 최근접(KNN) 그래프를 만들고, 
# 이를 바탕으로 Louvain 알고리즘을 사용해 세포그룹(cluster)을 나눈다.

pbmc <- FindNeighbors(pbmc, dims = 1:10) 
# - dims : PCA결과에서 몇 개의 주성분(PC)를 기반으로 이웃을 계산할지 지정
# - PCA1번부터 10번까지의 축만 사용해서 세포들 간 유사도를 계산하겠다는 뜻! 

# dims = 1:15 또는 dims = 1:20 등으로 바꿔서 실험해보는 것도 가능합니다.


# 그래프를 cluster로 나누는 과정 
# 유사한 세포끼리 모여있는 커뮤니티 찾기 
pbmc <- FindClusters(pbmc, resolution = 0.5)

# 클러스터의 결과 확인하기 
# 각 세포가 몇 번 클러스터에 속해있는지 알려준다. 
head(Idents(pbmc), 5)


# ----
### Run non-linear dimensional reduction (UMAP/tSNE)
# 데이터 시각화하기
# UMAP/tSNE -> 시각화 도구로는 ㄱㅊ, 그러나 해석의 근거로만 쓰면 X

# UMAP 계산
# - 비슷한 세포끼리는 가까이, 다른 샘플은 멀리 보이도록 해주는 방식 
pbmc <- RunUMAP(pbmc, dims = 1:10)

# 계산된 UMAP 좌표를 시각화 한다. 
DimPlot(pbmc, reduction = "umap", label = TRUE)

# 저장하기
saveRDS(pbmc, file = "./output/pbmc_tutorial.rds")

# ----
### Finding differentially expressed features (cluster biomarkers)
# 클러스터 마커 찾기, 시각화 

# 각 클러스터의 대표 유전자를 찾아주는 기능 제공한다.
# - 특정 클러스터 마커 찾기 , 여기서는 iden.1 = 2이므로 2번 유전자 마터를 찾는다.
# - 나머지  셀들을 비교해서, 2번 클러스터에 특이적으로 발현된 유전자를 찾아준다. 
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n=5)

# avg_log2FC ㅖ: 클러스터간 유전자 발현 차이
# pct.1, pct.2 -> 유전자가 해당 클러스터와 비교 그룹에서발현한 세포 비율 

# 클러스터 5번을  0,3와 비교해서 클러스터 5를 대표하는 유전자를 찾는다. 
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0,3))
head(cluster5.markers, n = 5)

# FindAllMarkers() : 모든 클러스터 마커 자동 검색 
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)

# 특정 기준 넘는 마커만 보기 
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# 다른 테스트도 가능 
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

# 마커를 시각화 해보자 
# - VlnPlot : 유전자 분포 
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

# - 유전자 별로 각 클러스터에서 발현할 확률 분포 
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

# FeaturePlot() : UMAP에 유전자 발현 위치 표시
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))

# DoHeatMap : 클러스터별 마커 히트맵
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC >1) %>%
  slice_head(n=10) %>%
  ungroup() -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

#----
### Assigning cell type identity to clusters 
# 분석에서 나온 클러스터(숫자)들을 이미 알려진 특정 세포 타입을 나타내는 표준 마커를 보고 어떤 세포인지 이름을 붙인다. 
# 클러스터 id와 세포 유형을 연결하는 과정 

new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")

# 클러스터 id와 이름 연결하기 
names(new.cluster.ids) <- levels(pbmc)

# cluster identity 이름 변경 
# - RenameIdents : 기존 클러스터 번호 대신 위에서 만든 세포 타입 이름을 적용한다. 
pbmc <- RenameIdents(pbmc, new.cluster.ids)

DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# 좀 더 꾸민 플롯 생성 및 저장

library(ggplot2)
plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
ggsave(filename = "./output/pbmc3k_umap.jpg", height = 7, width = 12, plot = plot, quality = 50)

saveRDS(pbmc, file = "./output/pbmc3k_final.rds")
