library(shiny)
library(ggplot2)

# 配置参数
#input_file <- "nodes.txt"    # 输入文件名
#output_file <- "output.pdf"  # 输出文件名
total_width <- 10            # 总宽度（单位：inch）
line_y <- 0                  # 主线Y坐标
text_offset <- -0.1            # 文字标注偏移量

# 绘图函数
myplot <- function(df) {
  # 按原始坐标排序
  sorted_df <- df[order(df$pos),]
  
  # 处理坐标标准化
  x_min <- min(sorted_df$pos)
  x_max <- max(sorted_df$pos)
  
  # 避免除零错误
  sorted_df$normalized <- ifelse(rep(x_max == x_min, dim(sorted_df)[1]), 
                                 1, 
                                 1 + (sorted_df$pos - x_min) * (total_width - 1) / (x_max - x_min))
  
  # 创建图形
  p <- ggplot(data=sorted_df) +
    # 绘制主线
    geom_segment(y=line_y, yend=line_y, x=1, xend=total_width, colour="black", size=1) +
    # 绘制节点标记
    geom_segment(aes(x=normalized, xend=normalized), y=line_y, yend=line_y + 0.1,
                 colour="black", size=1) +
    # 添加文字标注
    geom_text(aes(x=normalized, y=line_y + text_offset, label=node), 
              angle=45, hjust=1, vjust=1, size=3) +
    # 美化图形
    coord_cartesian(xlim=c(0, total_width + 1), ylim=c(line_y - 1, line_y + 0.5)) +
    theme_void() + # 移除所有轴和背景
    geom_blank()
  
  return(p)
}

# Define UI ----
ui <- fluidPage(
  titlePanel("GGplot2绘图应用"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "选择nodes.txt文件",
                accept = c("text/plain", ".txt")),
      textInput("outputFilename", "输出文件名", value = "output")
    ),
    mainPanel(
      plotOutput("plot"),
      downloadButton("downloadPlot", "下载图形")
    )
  )
)

# Define server logic ----
server <- function(input, output, session) {
  # 读取上传的文件
  data <- reactive({
    req(input$file1)
    read.table(input$file1$datapath, header = F, col.names=c('node', 'pos'))
  })
  
  # 生成并显示图形
  # 注意!output$plot仅用于预览图片，不能直接用来下载！
  output$plot <- renderPlot({
    myplot(data())
  })
  
  # 提供下载图形的链接
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste0(input$outputFilename, ".pdf")
    },
    content = function(file) {
      p <- myplot(data())
      ggsave(file, plot = p, device = "pdf",
             width=total_width+1, height=2, unit="in", dpi=300)
    }
  )
}

# Run the app ----
shinyApp(ui = ui, server = server)
