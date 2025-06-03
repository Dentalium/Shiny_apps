library(shiny)
library(ggplot2)
library(tidyr)
library(dplyr)

# 自定义函数：核心计算逻辑
calculate_qpcr <- function(df, ref.group, ref.gene, group_levels = NULL) {

  # 检查对照组是否存在
  if (!ref.group %in% df$group) {
    stop("错误：对照组 '", ref.group, "' 未在数据中找到！", head(df))
  }
  
  # 检查内参基因是否存在
  if (!ref.gene %in% df$gene) {
    stop("错误：内参基因 '", ref.gene, "' 未在数据中找到！")
  }
  # res
  res.all <- NULL
  
  # 目的基因
  target.genes <- dplyr::setdiff(unique(df$gene), ref.gene)
  
  #genes <- target.genes[1] # 只测试一个基因，分步运行
  for (genes in target.genes) {
    #{
    # 提取目标基因和内参基因
    df.sub <- df %>%
      dplyr::filter(gene %in% c(genes, ref.gene))
    
    # 提取对照组
    df.sub.ck <- df.sub %>%
      dplyr::filter(group == ref.group)
    
    # reference gene in CK
    # 提取对照组内参（所有组）
    df.sub.ck.ref.gene <- df.sub.ck %>%
      dplyr::filter(gene == ref.gene)
    # 计算对照组内参的平均值
    mean.ck.ref.gene <- mean(df.sub.ck.ref.gene$cq)
    
    # 提取对照组目标基因（所有组）
    df.sub.ck.target.gene <- df.sub.ck %>%
      dplyr::filter(gene != ref.gene)
    # 计算对照组目标基因平均值
    mean.ck.target.gene <- mean(df.sub.ck.target.gene$cq)
    
    # 目标基因-内参（对照组）
    dct0 <- mean.ck.target.gene - mean.ck.ref.gene
    
    # 遍历所有组（包括对照组）
    for (groups in unique(df.sub$group)) {
      # groups <- unique(df.sub$group)[1]  # 只测试一个组(对照组)
      #{
      # 提取当前组的目的基因和内参，技术重复之间取平均值
      df.sub.group <- df.sub %>%
        dplyr::filter(group == groups) %>%
        dplyr::select(biorep, gene, cq) %>%
        dplyr::group_by(gene, biorep) %>%
        dplyr::mutate(cq = mean(cq)) %>%
        dplyr::ungroup() %>%
        dplyr::distinct_all() %>%
        tidyr::pivot_wider(id_cols = "biorep", names_from = "gene", values_from = "cq") %>%
        dplyr::mutate(dct0 = dct0)
      # including ref gene or not
      # 如果提取出的数据不包含目的基因则报错
      if (ncol(df.sub.group) == 3 & !genes %in% colnames(df.sub.group)) {
        stop(paste0("Data of target gene ", genes, " has some problem, please check it and try again!"))
      }
      # 选择正确的内参
      if (colnames(df.sub.group)[3] == ref.gene) { # 如果第三列是内参
        new_name <- c("biorep", "Target", "Reference", "dct0")
      } else { # 如果不是
        new_name <- c("biorep", "Reference", "Target", "dct0")
      }
      df.sub.group %>%
        #magrittr::set_names(c("biorep", "Target", "Reference", "ddct1")) %>%
        dplyr::rename(!!!setNames(names(df.sub.group), new_name)) |>
        dplyr::mutate(expression = 2^(-(Target - Reference - dct0))) %>% # 每个生物学分别ddct计算表达量
        dplyr::mutate(
          group = groups,
          gene = genes
        ) %>%
        dplyr::select(group, gene, biorep, expression) %>%
        rbind(res.all) -> res.all
    }
  }
  
  # 统计汇总
  res_summary <- res.all %>%
    group_by(group, gene) %>%
    summarize(mean=mean(expression),
              sd=sd(expression))
  # 因子水平
  myfactors <- group_levels
  if (!is.null(myfactors)) {
    res_summary$group <- factor(res_summary$group, levels = myfactors)
    print(myfactors)
  }

  return(list(
    raw_result = res.all,
    summary_stats = res_summary
  ))
}

plot_qpcr <- function(df) {
  p <- ggplot(data = df) +
    geom_col(aes(x=group, y=mean, fill=gene), position="dodge", width=0.7) +
    # 误差棒
    geom_errorbar(aes(x=group, ymin=mean-sd, ymax=mean+sd, group=gene),
                  position=position_dodge(width = 0.7), width=0.25, linewidth=0.25) +
    scale_fill_viridis_d(begin = 0.2, end = 0.9, direction = -1) +
    labs(y="expression")
  
  return(p)
}

# UI界面
ui <- fluidPage(
  titlePanel("qPCR定量工具"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "上传原始数据 (TXT)", accept = ".txt"),
      textInput("control", "对照组名称", placeholder = "例如: CK"),
      textInput("housekeeping", "内参基因名称", placeholder = "例如: OsUBQ"),
      textInput("groups", "处理组顺序 (逗号分隔)(选填)", placeholder = "例如: Treatment,CK"),
      actionButton("analyze", "开始分析", class = "btn-primary"),
      tags$hr(),
      downloadButton("download_raw", "下载定量数据"),
      downloadButton("download_stats", "下载统计结果"),
      downloadButton("download_plot", "下载图表"),
      downloadLink("example", "下载示例数据"),
      helpText("原始数据要求："),
      helpText("1. 第一行为列名，列名依次是id、cq、group、gene、biorep；"),
      helpText("2. id为样本编号，可以任意填写；group为处理/分组；gene为基因名称，biorep为生物学重复；"),
      helpText("3. 文件必须是制表符分隔的文本文件，必须是UTF-8编码，请正确使用excel“另存为”功能。")
      
    ),
    mainPanel(
      plotOutput("expr_plot"),
      h4("定量数据预览"),
      tableOutput("raw_table"),
      h4("统计结果预览"),
      tableOutput("stats_table")
    )
  )
)

# 服务器逻辑
server <- function(input, output) {
  
  # 解析分组顺序
  group_levels <- reactive({
    if(nchar(input$groups) > 0) {
      unlist(strsplit(input$groups, ",")) %>% trimws()
    } else {
      NULL
    }
  })
  
  # 读取上传数据
  raw_data <- reactive({
    req(input$file)
    read.delim(input$file$datapath, stringsAsFactors = FALSE)
  })
  
  # 执行分析
  analysis_results <- eventReactive(input$analyze, {
    req(raw_data(), input$control, input$housekeeping)
    
    # 验证内参基因存在
    #validate(
    #  need(input$housekeeping %in% raw_data()$gene,
    #       "错误：输入的内参基因名称未在数据中找到！")
    #)
    
    calculate_qpcr(
      df = raw_data(),
      ref.group = input$control,
      ref.gene = input$housekeeping,
      group_levels = group_levels()
    )
  })
  
  # 绘制图表
  plot_object <- reactive({
    res_summ <- analysis_results()$summary_stats
    plot_qpcr(df = res_summ)
  })
  # 预览
  output$expr_plot <- renderPlot({
    res_summ <- analysis_results()$summary_stats
    plot_object()
  })
  
  # 数据预览
  output$raw_table <- renderTable({
    head(analysis_results()$raw_result)
  })
  
  output$stats_table <- renderTable({
    analysis_results()$summary_stats
  })
  
  # 下载处理
  output$download_raw <- downloadHandler(
    filename = function() {"qpcr_raw_results.csv"},
    content = function(file) {
      write.csv(analysis_results()$raw_result, file, row.names = FALSE)
    }
  )
  
  output$download_stats <- downloadHandler(
    filename = function() {"qpcr_summary_stats.csv"},
    content = function(file) {
      write.csv(analysis_results()$summary_stats, file, row.names = FALSE)
    }
  )
  
  output$download_plot <- downloadHandler(
    filename = function() {"expression_plot.pdf"},
    content = function(file) {
      ggsave(file, plot = plot_object(), device = "pdf", 
             width = 10, height = 6, dpi = 300)
    }
  )
  
  output$example <- downloadHandler(
    filename = function() {"example.txt"},
    content = function(file) {
      file_path <- normalizePath(file.path("example_data", "example.txt"))
      validate(need(file.exists(file_path), "示例文件未找到！"))
      file.copy(file_path, file)
    }
  )
}

shinyApp(ui = ui, server = server)
