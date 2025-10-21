# ============================================================================
# Server Logic for Drug Splicing Effects Explorer
# ============================================================================

server <- function(input, output, session) {
  
  # ============================================================================
  # REACTIVE DATA
  # ============================================================================
  
  # Get filtered datasets based on sidebar filters
  filtered_data <- reactive({
    req(datasets, input$fdr_threshold, input$deltapsi_threshold)
    
    # Map drug names to dataset keys
    drug_map <- c("Tubercidin" = "tub", "Pladienolide B" = "pladb", "Spliceostatin A" = "ssa")
    
    filtered <- list()
    for (drug in input$selected_drugs) {
      key <- drug_map[drug]
      if (!is.null(datasets[[key]])) {
        # Match original code: na.omit BEFORE filtering
        data <- datasets[[key]] %>%
          na.omit() %>%  # Remove rows with any NA values (like original code)
          mutate(
            event_type = classify_event(EVENT),
            significant = FDR <= input$fdr_threshold & abs(deltapsi) >= input$deltapsi_threshold,
            direction = ifelse(deltapsi > 0, "Increased", "Decreased")
          ) %>%
          filter(event_type %in% input$event_types)
        
        filtered[[drug]] <- data
      }
    }
    filtered
  })
  
  # Get significant events only
  significant_data <- reactive({
    req(filtered_data())
    lapply(filtered_data(), function(x) filter(x, significant == TRUE))
  })
  
  # ============================================================================
  # OVERVIEW TAB
  # ============================================================================
  
  # Value boxes for each drug
  output$vbox_tub <- renderValueBox({
    data <- filtered_data()[["Tubercidin"]]
    if (is.null(data)) {
      n_sig <- 0
      total <- 0
    } else {
      n_sig <- sum(data$significant, na.rm = TRUE)
      total <- nrow(data)
    }
    
    valueBox(
      value = format(n_sig, big.mark = ","),
      subtitle = paste("Significant Events - Tubercidin"),
      icon = icon("flask"),
      color = "green"
    )
  })
  
  output$vbox_pladb <- renderValueBox({
    data <- filtered_data()[["Pladienolide B"]]
    if (is.null(data)) {
      n_sig <- 0
      total <- 0
    } else {
      n_sig <- sum(data$significant, na.rm = TRUE)
      total <- nrow(data)
    }
    
    valueBox(
      value = format(n_sig, big.mark = ","),
      subtitle = paste("Significant Events - Pladienolide B"),
      icon = icon("flask"),
      color = "aqua"
    )
  })
  
  output$vbox_ssa <- renderValueBox({
    data <- filtered_data()[["Spliceostatin A"]]
    if (is.null(data)) {
      n_sig <- 0
      total <- 0
    } else {
      n_sig <- sum(data$significant, na.rm = TRUE)
      total <- nrow(data)
    }
    
    valueBox(
      value = format(n_sig, big.mark = ","),
      subtitle = paste("Significant Events - Spliceostatin A"),
      icon = icon("flask"),
      color = "purple"
    )
  })
  
  # Summary statistics table
  output$summary_table <- renderDT({
    req(filtered_data())
    
    summary_list <- lapply(names(filtered_data()), function(drug_name) {
      data <- filtered_data()[[drug_name]]
      sig_data <- filter(data, significant == TRUE)
      
      data.frame(
        Drug = drug_name,
        `Total Events` = nrow(data),
        `Significant Events` = nrow(sig_data),
        `% Significant` = round(100 * nrow(sig_data) / nrow(data), 2),
        `Increased Inclusion` = sum(sig_data$deltapsi > 0, na.rm = TRUE),
        `Decreased Inclusion` = sum(sig_data$deltapsi < 0, na.rm = TRUE),
        `Median |ΔPSI|` = round(median(abs(sig_data$deltapsi), na.rm = TRUE), 3),
        `Max |ΔPSI|` = round(max(abs(sig_data$deltapsi), na.rm = TRUE), 3),
        check.names = FALSE
      )
    })
    
    summary_df <- bind_rows(summary_list)
    
    datatable(
      summary_df,
      options = list(
        dom = 't',
        pageLength = 3,
        ordering = FALSE
      ),
      rownames = FALSE
    ) %>%
      formatStyle(
        'Drug',
        target = 'row',
        backgroundColor = styleEqual(
          c('Tubercidin', 'Pladienolide B', 'Spliceostatin A'),
          c('#A0C1B920', '#70A0AF20', '#70699320')
        )
      )
  })
  
  # ============================================================================
  # VOLCANO PLOT TAB
  # ============================================================================
  
  volcano_plot_data <- reactive({
    req(input$volcano_drug)
    data <- filtered_data()[[input$volcano_drug]]
    req(data)
    
    # Prepare data for plotting
    data <- data %>%
      mutate(
        negLog10FDR = pmin(-log10(pmax(FDR, 1e-300)), 4),  # Cap at 4 (FDR = 0.0001)
        color_group = case_when(
          significant & deltapsi > 0 ~ "Increased",
          significant & deltapsi < 0 ~ "Decreased",
          TRUE ~ "Not Significant"
        ),
        label = ifelse(significant, GENE, NA)
      )
    
    # Sample non-significant points for performance (keep all significant)
    sig_data <- data %>% filter(significant)
    nonsig_data <- data %>% filter(!significant)
    
    # If too many non-significant points, sample them
    if (nrow(nonsig_data) > 5000) {
      set.seed(42)
      nonsig_data <- sample_n(nonsig_data, 5000)
    }
    
    data <- bind_rows(sig_data, nonsig_data)
    
    # Get top events for labeling
    if (input$volcano_labels) {
      top_events <- sig_data %>%
        arrange(desc(abs(deltapsi))) %>%
        head(input$volcano_top_n)
      
      data$show_label <- data$EVENT %in% top_events$EVENT
    } else {
      data$show_label <- FALSE
    }
    
    data
  })
  
  output$volcano_plot <- renderPlotly({
    req(volcano_plot_data())
    data <- volcano_plot_data()
    
    # Color mapping - use actual color values, not named vector
    color_map <- c(
      "Increased" = unname(SPLICING_COLORS["Included"]),
      "Decreased" = unname(SPLICING_COLORS["Skipped"]),
      "Not Significant" = unname(SPLICING_COLORS["Unchanged"])
    )
    
    # Ensure color_group is a factor with correct levels
    data <- data %>%
      mutate(color_group = factor(color_group, levels = names(color_map)))
    
    # Base plot without labels (plotly doesn't support geom_text_repel)
    p <- ggplot(data, aes(x = deltapsi, y = negLog10FDR, 
                          color = color_group,
                          text = paste0(
                            "Gene: ", GENE, "\n",
                            "Event: ", EVENT, "\n",
                            "ΔPSI: ", round(deltapsi, 3), "\n",
                            "FDR: ", format.pval(FDR, digits = 2)
                          ))) +
      geom_point(aes(size = negLog10FDR), alpha = 0.6) +
      scale_color_manual(values = color_map, name = "", drop = FALSE) +
      scale_size_continuous(range = c(1, 4), guide = "none") +
      geom_hline(yintercept = -log10(input$fdr_threshold), 
                 linetype = "dashed", color = "gray50", linewidth = 0.5) +
      geom_vline(xintercept = c(-input$deltapsi_threshold, input$deltapsi_threshold), 
                 linetype = "dashed", color = "gray50", linewidth = 0.5) +
      labs(
        title = paste("Volcano Plot -", input$volcano_drug),
        x = "ΔPSI (Change in Percent Spliced In)",
        y = "-log10(FDR)"
      ) +
      theme_publication()
    
    # Convert to plotly with enhanced tooltip
    ply <- ggplotly(p, tooltip = "text") %>%
      layout(
        legend = list(orientation = "h", y = -0.15),
        hovermode = "closest"
      )
    
    # Add annotations for top genes if labels requested
    if (input$volcano_labels && any(data$show_label)) {
      label_data <- filter(data, show_label)
      
      # Add text annotations (plotly native way)
      for (i in 1:nrow(label_data)) {
        ply <- ply %>%
          add_annotations(
            x = label_data$deltapsi[i],
            y = label_data$negLog10FDR[i],
            text = label_data$GENE[i],
            showarrow = TRUE,
            arrowhead = 2,
            arrowsize = 0.5,
            arrowwidth = 1,
            arrowcolor = "gray30",
            ax = 20,
            ay = -30,
            font = list(size = 10, color = "black"),
            bgcolor = "rgba(255, 255, 255, 0.8)",
            bordercolor = "gray",
            borderwidth = 1
          )
      }
    }
    
    ply
  })
  
  # Volcano direction plot
  output$volcano_direction_plot <- renderPlot({
    req(volcano_plot_data())
    data <- volcano_plot_data() %>% filter(significant)
    
    if (nrow(data) == 0) {
      plot.new()
      text(0.5, 0.5, "No significant events", cex = 1.5)
      return()
    }
    
    direction_counts <- data %>%
      count(direction) %>%
      mutate(direction = factor(direction, levels = c("Increased", "Decreased")))
    
    # Use unname to avoid named vector issues
    fill_colors <- c(
      "Increased" = unname(SPLICING_COLORS["Included"]),
      "Decreased" = unname(SPLICING_COLORS["Skipped"])
    )
    
    ggplot(direction_counts, aes(x = direction, y = n, fill = direction)) +
      geom_col(width = 0.7) +
      geom_text(aes(label = n), vjust = -0.5, size = 5, fontface = "bold") +
      scale_fill_manual(values = fill_colors, drop = FALSE) +
      labs(
        title = "Significant Events by Direction",
        x = NULL,
        y = "Count"
      ) +
      theme_publication() +
      theme(legend.position = "none") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
  })
  
  # Event types distribution in volcano
  output$volcano_event_types_plot <- renderPlot({
    req(volcano_plot_data())
    data <- volcano_plot_data() %>% filter(significant)
    
    if (nrow(data) == 0) {
      plot.new()
      text(0.5, 0.5, "No significant events", cex = 1.5)
      return()
    }
    
    event_counts <- data %>%
      count(event_type) %>%
      arrange(desc(n))
    
    # Create color mapping for only the event types present in data
    event_colors_present <- EVENT_COLORS[names(EVENT_COLORS) %in% event_counts$event_type]
    event_colors_present <- unname(event_colors_present)
    names(event_colors_present) <- names(EVENT_COLORS)[names(EVENT_COLORS) %in% event_counts$event_type]
    
    ggplot(event_counts, aes(x = reorder(event_type, n), y = n, fill = event_type)) +
      geom_col(width = 0.7) +
      geom_text(aes(label = n), hjust = -0.2, size = 4, fontface = "bold") +
      scale_fill_manual(values = event_colors_present, drop = TRUE) +
      coord_flip() +
      labs(
        title = "Event Type Distribution",
        x = NULL,
        y = "Count"
      ) +
      theme_publication() +
      theme(legend.position = "none") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
  })
  
  # ============================================================================
  # EVENT ANALYSIS TAB
  # ============================================================================
  
  output$events_barplot <- renderPlot({
    req(significant_data())
    
    # Prepare data
    event_data <- lapply(names(significant_data()), function(drug) {
      data <- significant_data()[[drug]]
      data %>%
        mutate(Drug = drug) %>%
        count(Drug, event_type)
    }) %>% bind_rows()
    
    if (nrow(event_data) == 0) {
      plot.new()
      text(0.5, 0.5, "No significant events", cex = 1.5)
      return()
    }
    
    # Order event types
    event_data <- event_data %>%
      mutate(event_type = factor(event_type, levels = c("Exon", "Intron", "Alt5", "Alt3")))
    
    # Create color mapping for only the event types present in data
    event_colors_present <- unname(EVENT_COLORS[names(EVENT_COLORS) %in% unique(event_data$event_type)])
    names(event_colors_present) <- names(EVENT_COLORS)[names(EVENT_COLORS) %in% unique(event_data$event_type)]
    
    if (input$event_view_type == "proportion") {
      event_data <- event_data %>%
        group_by(Drug) %>%
        mutate(prop = n / sum(n) * 100) %>%
        ungroup()
      
      p <- ggplot(event_data, aes(x = Drug, y = prop, fill = event_type)) +
        geom_col(position = "stack", width = 0.7) +
        labs(y = "Percentage (%)", title = "Splicing Event Distribution (Proportions)")
      
      if (input$event_show_numbers) {
        p <- p + geom_text(
          aes(label = sprintf("%.1f%%", prop)),
          position = position_stack(vjust = 0.5),
          size = 4,
          color = "white",
          fontface = "bold"
        )
      }
    } else {
      p <- ggplot(event_data, aes(x = Drug, y = n, fill = event_type)) +
        geom_col(position = "dodge", width = 0.7) +
        labs(y = "Count", title = "Splicing Event Counts")
      
      if (input$event_show_numbers) {
        p <- p + geom_text(
          aes(label = n),
          position = position_dodge(width = 0.7),
          vjust = -0.5,
          size = 4,
          fontface = "bold"
        )
      }
    }
    
    p + 
      scale_fill_manual(values = event_colors_present, name = "Event Type", drop = FALSE) +
      theme_publication() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right"
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
  })
  
  output$events_summary_table <- renderDT({
    req(significant_data())
    
    summary_list <- lapply(names(significant_data()), function(drug) {
      data <- significant_data()[[drug]]
      
      event_summary <- data %>%
        count(event_type) %>%
        pivot_wider(names_from = event_type, values_from = n, values_fill = 0) %>%
        mutate(Drug = drug, Total = rowSums(across(where(is.numeric)))) %>%
        select(Drug, Total, everything())
    }) %>% bind_rows()
    
    datatable(
      summary_list,
      options = list(dom = 't', ordering = FALSE),
      rownames = FALSE
    )
  })
  
  output$events_pie_chart <- renderPlot({
    req(input$events_pie_drug, significant_data())
    
    # Check if drug exists in significant_data
    if (!input$events_pie_drug %in% names(significant_data())) {
      plot.new()
      text(0.5, 0.5, "No data available for selected drug", cex = 1.5)
      return()
    }
    
    data <- significant_data()[[input$events_pie_drug]]
    
    if (is.null(data) || nrow(data) == 0) {
      plot.new()
      text(0.5, 0.5, "No significant events", cex = 1.5)
      return()
    }
    
    pie_data <- data %>%
      count(event_type) %>%
      mutate(
        percentage = n / sum(n) * 100,
        label = paste0(event_type, "\n", n, " (", round(percentage, 1), "%)")
      )
    
    # Create color mapping for only the event types present
    event_colors_present <- unname(EVENT_COLORS[names(EVENT_COLORS) %in% pie_data$event_type])
    names(event_colors_present) <- names(EVENT_COLORS)[names(EVENT_COLORS) %in% pie_data$event_type]
    
    ggplot(pie_data, aes(x = "", y = n, fill = event_type)) +
      geom_col(width = 1, color = "white", linewidth = 1) +
      coord_polar(theta = "y") +
      scale_fill_manual(values = event_colors_present, drop = TRUE) +
      labs(title = paste("Event Distribution -", input$events_pie_drug)) +
      theme_void() +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        legend.position = "right",
        legend.title = element_blank()
      ) +
      geom_text(
        aes(label = paste0(round(percentage, 1), "%")),
        position = position_stack(vjust = 0.5),
        size = 5,
        fontface = "bold",
        color = "white"
      )
  })
  
  # ============================================================================
  # UPSET PLOT TAB
  # ============================================================================
  
  upset_data <- reactive({
    req(significant_data())
    
    if (input$upset_level == "genes") {
      gene_lists <- lapply(names(significant_data()), function(drug_name) {
        data <- significant_data()[[drug_name]]
        if (input$upset_event_filter != "all") {
          data <- filter(data, event_type == input$upset_event_filter)
        }
        unique(data$GENE)
      })
      names(gene_lists) <- names(significant_data())
    } else {
      gene_lists <- lapply(names(significant_data()), function(drug_name) {
        data <- significant_data()[[drug_name]]
        if (input$upset_event_filter != "all") {
          data <- filter(data, event_type == input$upset_event_filter)
        }
        unique(data$EVENT)
      })
      names(gene_lists) <- names(significant_data())
    }
    
    gene_lists
  })
  
  output$upset_plot <- renderPlot({
    req(upset_data())
    
    sets <- upset_data()
    if (length(sets) < 2) {
      plot.new()
      text(0.5, 0.5, "Need at least 2 drugs selected", cex = 1.5)
      return()
    }
    
    # Get all unique items
    all_items <- unique(unlist(sets))
    
    # Create combinations with meaningful labels
    n_sets <- length(sets)
    set_names <- names(sets)
    
    # Calculate all possible intersections
    intersections <- list()
    
    # Single sets (unique to each)
    for (i in 1:n_sets) {
      unique_items <- sets[[i]]
      for (j in setdiff(1:n_sets, i)) {
        unique_items <- setdiff(unique_items, sets[[j]])
      }
      if (length(unique_items) > 0) {
        intersections[[paste0("Only ", set_names[i])]] <- length(unique_items)
      }
    }
    
    # Pairwise intersections (but not in third if applicable)
    if (n_sets == 2) {
      shared <- intersect(sets[[1]], sets[[2]])
      if (length(shared) > 0) {
        intersections[[paste(set_names[1], "∩", set_names[2])]] <- length(shared)
      }
    } else if (n_sets == 3) {
      # Two-way intersections (excluding third)
      for (i in 1:(n_sets-1)) {
        for (j in (i+1):n_sets) {
          shared <- intersect(sets[[i]], sets[[j]])
          third <- setdiff(1:n_sets, c(i, j))
          shared <- setdiff(shared, sets[[third]])
          if (length(shared) > 0) {
            intersections[[paste(set_names[i], "∩", set_names[j])]] <- length(shared)
          }
        }
      }
      
      # Three-way intersection
      shared_all <- Reduce(intersect, sets)
      if (length(shared_all) > 0) {
        intersections[[paste(set_names, collapse = " ∩ ")]] <- length(shared_all)
      }
    }
    
    # Convert to data frame
    if (length(intersections) == 0) {
      plot.new()
      text(0.5, 0.5, "No overlaps found", cex = 1.5)
      return()
    }
    
    plot_data <- data.frame(
      category = names(intersections),
      count = unlist(intersections)
    ) %>%
      arrange(desc(count))
    
    # Create bar plot
    ggplot(plot_data, aes(x = reorder(category, count), y = count)) +
      geom_col(fill = "#706993", width = 0.7) +
      geom_text(aes(label = count), hjust = -0.2, size = 5, fontface = "bold") +
      coord_flip() +
      labs(
        title = "Gene/Event Set Overlaps",
        x = NULL,
        y = "Count"
      ) +
      theme_publication() +
      theme(axis.text.y = element_text(size = 11)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
  })
  
  output$upset_stats_table <- renderDT({
    req(upset_data())
    
    sets <- upset_data()
    set_names <- names(sets)
    
    # Calculate statistics
    stats_list <- lapply(set_names, function(name) {
      total_items <- length(sets[[name]])
      
      # Calculate unique (not in any other set)
      unique_items <- sets[[name]]
      for (other in setdiff(set_names, name)) {
        unique_items <- setdiff(unique_items, sets[[other]])
      }
      
      data.frame(
        Drug = name,
        `Total Items` = total_items,
        `Unique Items` = length(unique_items),
        check.names = FALSE
      )
    })
    
    stats_df <- bind_rows(stats_list)
    
    # Add pairwise overlaps as additional columns
    if (length(set_names) >= 2) {
      for (i in 1:(length(set_names)-1)) {
        for (j in (i+1):length(set_names)) {
          overlap <- length(intersect(sets[[set_names[i]]], sets[[set_names[j]]]))
          col_name <- paste0(set_names[i], " ∩ ", set_names[j])
          stats_df[[col_name]] <- c(overlap, rep(NA, nrow(stats_df) - 1))
        }
      }
    }
    
    # Add three-way overlap if applicable
    if (length(set_names) == 3) {
      all_overlap <- length(Reduce(intersect, sets))
      stats_df[["All Three"]] <- c(all_overlap, rep(NA, nrow(stats_df) - 1))
    }
    
    datatable(
      stats_df,
      options = list(dom = 't', ordering = FALSE),
      rownames = FALSE
    )
  })
  
  output$venn_plot <- renderPlot({
    req(upset_data())
    
    sets <- upset_data()
    
    if (length(sets) < 2 || length(sets) > 3) {
      plot.new()
      text(0.5, 0.5, "Venn diagram available for 2-3 drugs only", cex = 1.5)
      return()
    }
    
    # Simple Venn representation using base graphics
    plot.new()
    plot.window(xlim = c(0, 10), ylim = c(0, 10))
    
    set_names <- names(sets)
    
    if (length(sets) == 2) {
      # Two circles
      symbols(c(3.5, 6.5), c(5, 5), circles = c(2, 2), 
              inches = FALSE, add = TRUE, fg = "black", lwd = 2)
      
      # Labels
      text(2, 5, set_names[1], cex = 1.2, font = 2)
      text(8, 5, set_names[2], cex = 1.2, font = 2)
      
      # Counts
      only_1 <- length(setdiff(sets[[1]], sets[[2]]))
      only_2 <- length(setdiff(sets[[2]], sets[[1]]))
      both <- length(intersect(sets[[1]], sets[[2]]))
      
      text(2.5, 5, only_1, cex = 1.5, font = 2)
      text(5, 5, both, cex = 1.5, font = 2)
      text(7.5, 5, only_2, cex = 1.5, font = 2)
      
    } else if (length(sets) == 3) {
      # Three circles
      symbols(c(3.5, 6.5, 5), c(6, 6, 3.5), circles = c(2, 2, 2), 
              inches = FALSE, add = TRUE, fg = "black", lwd = 2)
      
      # Labels
      text(2, 7.5, set_names[1], cex = 1, font = 2)
      text(8, 7.5, set_names[2], cex = 1, font = 2)
      text(5, 1.5, set_names[3], cex = 1, font = 2)
      
      # Calculate all regions
      all_three <- Reduce(intersect, sets)
      only_1_2 <- setdiff(intersect(sets[[1]], sets[[2]]), sets[[3]])
      only_1_3 <- setdiff(intersect(sets[[1]], sets[[3]]), sets[[2]])
      only_2_3 <- setdiff(intersect(sets[[2]], sets[[3]]), sets[[1]])
      only_1 <- setdiff(setdiff(sets[[1]], sets[[2]]), sets[[3]])
      only_2 <- setdiff(setdiff(sets[[2]], sets[[1]]), sets[[3]])
      only_3 <- setdiff(setdiff(sets[[3]], sets[[1]]), sets[[2]])
      
      # Place counts
      text(2.5, 6, length(only_1), cex = 1.2, font = 2)
      text(7.5, 6, length(only_2), cex = 1.2, font = 2)
      text(5, 2.5, length(only_3), cex = 1.2, font = 2)
      text(5, 6, length(only_1_2), cex = 1.2, font = 2)
      text(3.5, 4, length(only_1_3), cex = 1.2, font = 2)
      text(6.5, 4, length(only_2_3), cex = 1.2, font = 2)
      text(5, 4.5, length(all_three), cex = 1.2, font = 2)
    }
    
    title(main = "Venn Diagram", cex.main = 1.5, font.main = 2)
  })
  
  # ============================================================================
  # HEATMAP TAB
  # ============================================================================
  
  output$heatmap_plot <- renderPlotly({
    req(input$heatmap_drug, input$heatmap_event_type)
    
    # Map drug name to key
    drug_map <- c("Tubercidin" = "tub", "Pladienolide B" = "pladb", "Spliceostatin A" = "ssa")
    drug_key <- drug_map[input$heatmap_drug]
    
    # Get significant events for the selected drug only
    data <- significant_data()[[input$heatmap_drug]]
    req(data)
    
    # Filter by event type
    if (input$heatmap_event_type == "exon") {
      data <- filter(data, grepl("EX", EVENT))
    } else {
      data <- filter(data, grepl("INT", EVENT))
    }
    
    # Limit to max events
    data <- data %>%
      arrange(desc(abs(deltapsi))) %>%
      head(input$heatmap_max_events)
    
    if (nrow(data) == 0) {
      return(plot_ly() %>% add_text(text = "No events match the selected criteria", x = 0.5, y = 0.5))
    }
    
    # Get PSI columns: control + selected drug only
    psi_cols <- grep("^reads_", names(data), value = TRUE)
    
    # Keep control columns (ctl) and selected drug columns
    drug_pattern <- paste0("_", tolower(drug_key), "_")
    psi_cols_filtered <- psi_cols[grepl("_ctl_", psi_cols, ignore.case = TRUE) | 
                          grepl(drug_pattern, psi_cols, ignore.case = TRUE)]
    
    # Sort columns: control first, then drug (by numeric prefix to maintain replicate order)
    psi_cols_ctl <- psi_cols_filtered[grepl("_ctl_", psi_cols_filtered)]
    psi_cols_drug <- psi_cols_filtered[grepl(drug_pattern, psi_cols_filtered)]
    
    # Sort each group by the numeric prefix (reads_N_...)
    psi_cols_ctl <- psi_cols_ctl[order(as.numeric(gsub("^reads_(\\d+)_.*", "\\1", psi_cols_ctl)))]
    psi_cols_drug <- psi_cols_drug[order(as.numeric(gsub("^reads_(\\d+)_.*", "\\1", psi_cols_drug)))]
    
    # Combine: control first, then drug
    psi_cols <- c(psi_cols_ctl, psi_cols_drug)
    
    if (length(psi_cols) == 0) {
      return(plot_ly() %>% add_text(text = "PSI data not available", x = 0.5, y = 0.5))
    }
    
    # Prepare matrix for clustering
    psi_matrix_wide <- data %>%
      select(EVENT, GENE, all_of(psi_cols))
    
    # Create labels
    row_labels <- paste0(psi_matrix_wide$GENE, " (", psi_matrix_wide$EVENT, ")")
    
    # Extract numeric matrix
    psi_numeric <- as.matrix(psi_matrix_wide[, psi_cols])
    rownames(psi_numeric) <- row_labels
    
    # Remove rows with NAs
    valid_rows <- complete.cases(psi_numeric)
    if (sum(valid_rows) < 2) {
      return(plot_ly() %>% add_text(text = "Not enough complete data for clustering", x = 0.5, y = 0.5))
    }
    
    psi_numeric <- psi_numeric[valid_rows, ]
    
    # Perform hierarchical clustering on rows
    row_dist <- dist(psi_numeric, method = "euclidean")
    row_hclust <- hclust(row_dist, method = input$heatmap_cluster_method)
    row_order <- row_hclust$order
    
    # Reorder data
    psi_matrix_wide_ordered <- psi_matrix_wide[valid_rows, ][row_order, ]
    
    # Prepare column names - clean them properly with explicit sample type
    col_names_clean <- sapply(psi_cols, function(col_name) {
      # Check if it's control or drug
      if (grepl("_ctl_", col_name, ignore.case = TRUE)) {
        # Control sample
        replicate_num <- gsub("^reads_(\\d+)_.*", "\\1", col_name)
        return(paste0("Control ", replicate_num))
      } else {
        # Drug sample - extract drug name and details
        clean <- gsub("^reads_\\d+_", "", col_name)
        clean <- gsub("_\\d{2}_\\d{2}_\\d{2}$", "", clean)
        clean <- gsub("_", " ", clean)
        clean <- tools::toTitleCase(clean)
        return(clean)
      }
    })
    
    # Make names unique (in case of duplicates)
    col_names_clean <- make.unique(col_names_clean, sep = " ")
    
    # Extract matrix and labels
    plot_matrix <- as.matrix(psi_matrix_wide_ordered[, psi_cols])
    row_labels_ordered <- paste0(psi_matrix_wide_ordered$GENE, " (", psi_matrix_wide_ordered$EVENT, ")")
    
    # Create hover text matrix with properly cleaned names
    hover_text <- matrix("", nrow = nrow(plot_matrix), ncol = ncol(plot_matrix))
    for (i in 1:nrow(plot_matrix)) {
      for (j in 1:ncol(plot_matrix)) {
        hover_text[i, j] <- paste0(
          "Event: ", row_labels_ordered[i], "<br>",
          "Sample: ", col_names_clean[j], "<br>",
          "PSI: ", round(plot_matrix[i, j], 1), "%"
        )
      }
    }
    
    # Create plotly heatmap with splicing colors (Removed=blue/low PSI, Retained=red/high PSI)
    plot_ly(
      x = col_names_clean,
      y = row_labels_ordered,
      z = plot_matrix,
      type = "heatmap",
      colorscale = list(
        c(0, "#B85450"),      # Removed (low PSI, skipped) - muted coral
        c(0.5, "#EEEEEE"),    # Middle - light gray
        c(1, "#D4A574")       # Retained (high PSI, included) - warm sand
      ),
      zmin = 0,
      zmax = 100,
      hovertemplate = "%{text}<extra></extra>",
      text = hover_text,
      colorbar = list(title = "PSI (%)")
    ) %>%
      layout(
        title = paste("PSI Heatmap -", input$heatmap_drug, "vs Control (clustered by", input$heatmap_cluster_method, ")"),
        xaxis = list(title = "", tickfont = list(size = 10)),
        yaxis = list(title = "", tickfont = list(size = 8)),
        margin = list(l = 200, b = 100)
      )
  })
  
  # ============================================================================
  # GO ENRICHMENT TAB
  # ============================================================================
  
  go_filtered_data <- reactive({
    req(input$go_drug, input$go_ontology)
    
    drug_map <- c("Tubercidin" = "tub", "Pladienolide B" = "pladb", "Spliceostatin A" = "ssa")
    drug_key <- drug_map[input$go_drug]
    
    go_result <- go_data[[drug_key]][[input$go_ontology]]
    
    if (is.null(go_result) || nrow(go_result) == 0) {
      return(NULL)
    }
    
    # Filter and prepare
    go_result %>%
      filter(!is.na(p.adjust), p.adjust <= input$go_pval_cutoff) %>%
      arrange(p.adjust) %>%
      head(input$go_top_n)
  })
  
  output$go_plot <- renderPlot({
    data <- go_filtered_data()
    
    if (is.null(data) || nrow(data) == 0) {
      plot.new()
      text(0.5, 0.5, "No GO terms meet the significance threshold", cex = 1.5)
      return()
    }
    
    # Prepare data
    plot_data <- data %>%
      mutate(
        negLogP = -log10(p.adjust),
        Description = ifelse(
          nchar(as.character(Description)) > 50,
          paste0(substr(Description, 1, 47), "..."),
          as.character(Description)
        )
      ) %>%
      arrange(negLogP)
    
    # Determine color column
    if ("Count" %in% names(plot_data)) {
      color_var <- "Count"
    } else if ("count" %in% names(plot_data)) {
      color_var <- "count"
    } else {
      plot_data$count <- 1
      color_var <- "count"
    }
    
    ggplot(plot_data, aes(x = negLogP, y = reorder(Description, negLogP))) +
      geom_col(aes_string(fill = color_var), width = 0.7) +
      scale_fill_gradient(low = "#E8E3F5", high = "#706993", name = "Gene\nCount") +
      labs(
        title = paste("GO Enrichment -", input$go_drug),
        subtitle = c(BP = "Biological Process", MF = "Molecular Function", 
                     CC = "Cellular Component")[input$go_ontology],
        x = "-log10(Adjusted P-value)",
        y = NULL
      ) +
      theme_publication() +
      theme(
        axis.text.y = element_text(size = 10),
        legend.position = "right"
      ) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.1)))
  })
  
  output$go_table <- renderDT({
    data <- go_filtered_data()
    
    if (is.null(data) || nrow(data) == 0) {
      return(datatable(data.frame(Message = "No significant GO terms found")))
    }
    
    # Select relevant columns
    display_cols <- intersect(
      c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "Count"),
      names(data)
    )
    
    display_data <- data[, display_cols]
    
    datatable(
      display_data,
      options = list(
        pageLength = 15,
        scrollX = TRUE,
        dom = 'Bfrtip'
      ),
      rownames = FALSE
    ) %>%
      formatSignif(columns = c("pvalue", "p.adjust"), digits = 3)
  })
  
  # ============================================================================
  # SEQUENCE FEATURES TAB
  # ============================================================================
  
  gc_plot_data <- reactive({
    req(input$gc_event_type, event_info)
    
    # Get significant events from all selected drugs with direction info
    sig_data <- lapply(names(significant_data()), function(drug) {
      data <- significant_data()[[drug]]
      
      if (input$gc_event_type == "exon") {
        data <- filter(data, grepl("EX", EVENT))
      } else {
        data <- filter(data, grepl("INT", EVENT))
      }
      
      data %>%
        mutate(
          Drug = drug,
          splicing = if_else(deltapsi < 0, "Removed", "Retained")
        ) %>%
        select(EVENT, GENE, Drug, deltapsi, splicing)
    }) %>% bind_rows()
    
    # Get ALL significant events across all drugs
    all_sig_events <- unique(sig_data$EVENT)
    
    # Get unchanged pool: events NOT significant in ANY drug (shared pool)
    # Use the first drug's FDR file to get all events, then filter
    first_drug <- names(filtered_data())[1]
    all_events_data <- filtered_data()[[first_drug]]
    
    if (input$gc_event_type == "exon") {
      all_events_data <- filter(all_events_data, grepl("EX", EVENT))
    } else {
      all_events_data <- filter(all_events_data, grepl("INT", EVENT))
    }
    
    # Filter to events NOT significant in ANY drug
    unchanged_pool <- all_events_data %>%
      filter(!EVENT %in% all_sig_events)
    
    # Sample max 1000 from the shared pool
    if (nrow(unchanged_pool) > 1000) {
      unchanged_pool <- unchanged_pool %>% sample_n(1000)
    }
    
    # Replicate unchanged pool for each drug (so it appears in each facet)
    unchanged_data_list <- lapply(names(significant_data()), function(drug) {
      unchanged_pool %>%
        mutate(
          Drug = drug,  # Assign to same drug for faceting
          splicing = "Unchanged"
        ) %>%
        select(EVENT, GENE, Drug, splicing)
    })
    
    unchanged_data <- bind_rows(unchanged_data_list)
    
    # Combine
    all_data_combined <- bind_rows(sig_data, unchanged_data)
    
    # Join with event info to get Seq_A (for introns/exons) and LE_o (length)
    if (!is.null(event_info)) {
      all_data_combined <- all_data_combined %>%
        left_join(event_info %>% select(EVENT, Seq_A, LE_o), by = "EVENT") %>%
        filter(!is.na(Seq_A))
      
      # Calculate GC content from Seq_A
      all_data_combined$GC <- sapply(all_data_combined$Seq_A, function(seq) {
        if (is.na(seq)) return(NA)
        seq <- toupper(as.character(seq))
        chars <- strsplit(seq, "")[[1]]
        gc <- sum(chars %in% c("G", "C"))
        n <- length(chars)
        if (n == 0) return(NA)
        100 * gc / n
      })
      
      # Set factor levels
      all_data_combined <- all_data_combined %>%
        mutate(
          splicing = factor(splicing, levels = c("Retained", "Removed", "Unchanged")),
          Drug = factor(Drug, levels = names(DRUG_COLORS))
        )
    }
    
    all_data_combined
  })
  
  output$gc_plot <- renderPlot({
    data <- gc_plot_data()
    req(data, nrow(data) > 0)
    
    # Color palette matching original
    splicing_colors <- c(
      "Retained" = "#D4A574",
      "Removed" = "#B85450",
      "Unchanged" = "#CCCCCC"
    )
    
    # Define comparisons for statistical tests (like original)
    comparisons <- list(
      c("Retained", "Removed"),
      c("Retained", "Unchanged"),
      c("Removed", "Unchanged")
    )
    
    # Create plot with splicing on x-axis
    p <- ggplot(data, aes(x = splicing, y = GC, fill = splicing)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6, linewidth = 0.5) +
      geom_jitter(width = 0.15, alpha = 0.5, size = 1.2) +
      stat_compare_means(
        method = "wilcox.test",
        comparisons = comparisons,
        label = "p.signif",
        label.size = 4
      ) +
      scale_fill_manual(values = splicing_colors) +
      facet_wrap(~ Drug, nrow = 1) +
      labs(
        title = paste("GC Content Distribution -", 
                      ifelse(input$gc_event_type == "exon", "Exons", "Introns")),
        x = "Splicing Outcome",
        y = "GC Content (%)"
      ) +
      theme_publication() +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = rel(0.85)),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = rel(0.95)),
        panel.spacing = unit(1, "lines")
      )
    
    if (input$gc_show_stats) {
      # Add median labels
      stat_data <- data %>%
        group_by(Drug, splicing) %>%
        summarise(
          median = median(GC, na.rm = TRUE),
          .groups = "drop"
        )
      
      p <- p + geom_text(
        data = stat_data,
        aes(label = sprintf("%.1f", median), y = median),
        vjust = -0.5,
        size = 3,
        fontface = "bold"
      )
    }
    
    p
  })
  
  # Length distribution
  output$length_plot <- renderPlot({
    req(input$length_event_type)
    
    data <- gc_plot_data()  # Reuse same data which has LE_o
    req(data, nrow(data) > 0)
    
    # Filter by event type if needed
    if (input$length_event_type != input$gc_event_type) {
      # Need to regenerate data for different event type
      if (input$length_event_type == "exon") {
        event_filter <- "EX"
      } else {
        event_filter <- "INT"
      }
      
      # Regenerate for length plot with shared unchanged pool
      sig_data <- lapply(names(significant_data()), function(drug) {
        df <- significant_data()[[drug]]
        df <- filter(df, grepl(event_filter, EVENT))
        df %>%
          mutate(
            Drug = drug,
            splicing = if_else(deltapsi < 0, "Removed", "Retained")
          ) %>%
          select(EVENT, GENE, Drug, deltapsi, splicing)
      }) %>% bind_rows()
      
      # Get ALL significant events across all drugs
      all_sig_events <- unique(sig_data$EVENT)
      
      # Get unchanged pool: events NOT significant in ANY drug (shared pool)
      first_drug <- names(filtered_data())[1]
      all_events_data <- filtered_data()[[first_drug]]
      all_events_data <- filter(all_events_data, grepl(event_filter, EVENT))
      
      # Filter to events NOT significant in ANY drug
      unchanged_pool <- all_events_data %>%
        filter(!EVENT %in% all_sig_events)
      
      # Sample max 1000 from the shared pool
      if (nrow(unchanged_pool) > 1000) {
        unchanged_pool <- unchanged_pool %>% sample_n(1000)
      }
      
      # Replicate unchanged pool for each drug (so it appears in each facet)
      unchanged_list <- lapply(names(significant_data()), function(drug) {
        unchanged_pool %>%
          mutate(
            Drug = drug,  # Assign to same drug for faceting
            splicing = "Unchanged"
          ) %>%
          select(EVENT, GENE, Drug, splicing)
      })
      
      data <- bind_rows(sig_data, bind_rows(unchanged_list)) %>%
        left_join(event_info %>% select(EVENT, Seq_A, LE_o), by = "EVENT") %>%
        filter(!is.na(Seq_A)) %>%
        mutate(
          splicing = factor(splicing, levels = c("Retained", "Removed", "Unchanged")),
          Drug = factor(Drug, levels = names(DRUG_COLORS))
        )
    }
    
    # Color palette
    splicing_colors <- c(
      "Retained" = "#D4A574",
      "Removed" = "#B85450",
      "Unchanged" = "#CCCCCC"
    )
    
    # Define comparisons for statistical tests (like original)
    comparisons <- list(
      c("Retained", "Removed"),
      c("Retained", "Unchanged"),
      c("Removed", "Unchanged")
    )
    
    # Clip at 95th percentile for display (like original)
    threshold_95 <- quantile(data$LE_o, 0.95, na.rm = TRUE)
    data <- data %>%
      mutate(
        LE_o_display = pmin(LE_o, threshold_95),
        is_outlier = LE_o > threshold_95
      )
    
    p <- ggplot(data, aes(x = splicing, y = LE_o_display, fill = splicing)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6, linewidth = 0.5) +
      geom_jitter(width = 0.15, alpha = 0.5, size = 1.2) +
      stat_compare_means(
        method = "wilcox.test",
        comparisons = comparisons,
        label = "p.signif",
        label.size = 4
      ) +
      scale_fill_manual(values = splicing_colors) +
      facet_wrap(~ Drug, nrow = 1) +
      labs(
        title = paste("Length Distribution -", 
                      ifelse(input$length_event_type == "exon", "Exons", "Introns")),
        x = "Splicing Outcome",
        y = "Length (nt)"
      ) +
      theme_publication() +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = rel(0.85)),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = rel(0.95)),
        panel.spacing = unit(1, "lines")
      )
    
    # Add caption about clipping
    n_outliers <- sum(data$is_outlier, na.rm = TRUE)
    if (n_outliers > 0) {
      p <- p + labs(caption = paste0("Note: ", n_outliers, " outliers > ", 
                                      round(threshold_95), " nt clipped for display; statistics on full data"))
    }
    
    if (input$length_show_stats) {
      stat_data <- data %>%
        group_by(Drug, splicing) %>%
        summarise(
          median = median(LE_o, na.rm = TRUE),
          .groups = "drop"
        )
      
      p <- p + geom_text(
        data = stat_data,
        aes(label = comma(round(median)), y = pmin(median, threshold_95)),
        vjust = -0.5,
        size = 3,
        fontface = "bold"
      )
    }
    
    p
  })
  
  # ============================================================================
  # RBP ANALYSIS TAB
  # ============================================================================
  
  rbp_data <- reactive({
    req(input$rbp_drug, input$rbp_event_type, input$rbp_direction)
    
    # Map drug name to file prefix
    drug_map <- c("Tubercidin" = "tub", "Pladienolide B" = "pladb", "Spliceostatin A" = "ssa")
    drug_key <- drug_map[input$rbp_drug]
    
    # Map event type to file suffix
    event_map <- c("ex" = "a3", "int" = "int")  # Simplified mapping
    event_key <- event_map[input$rbp_event_type]
    
    # Determine file path based on direction
    if (input$rbp_direction == "up") {
      file_name <- "pVal.up.vs.bg.RNAmap.txt"
    } else if (input$rbp_direction == "down") {
      file_name <- "pVal.dn.vs.bg.RNAmap.txt"
    } else {
      file_name <- "pVal.up.vs.bg.RNAmap.txt"  # Default to up
    }
    
    # Construct path
    rbp_dir <- sprintf("../rmaps_out/%s_%s_rmaps", drug_key, event_key)
    rbp_file <- file.path(rbp_dir, file_name)
    
    # Try to load
    if (file.exists(rbp_file)) {
      tryCatch({
        data <- read.delim(rbp_file, header = TRUE, stringsAsFactors = FALSE)
        return(data)
      }, error = function(e) {
        return(NULL)
      })
    } else {
      return(NULL)
    }
  })
  
  output$rbp_plot <- renderPlot({
    data <- rbp_data()
    
    if (is.null(data)) {
      plot.new()
      text(0.5, 0.5, "RBP data not available for this selection.\nRNA map files not found.", cex = 1.2)
      return()
    }
    
    # Extract RBP names (remove motif pattern)
    data$RBP_clean <- gsub("\\..*", "", data$RBP)
    
    # Calculate average significance across regions
    pval_cols <- grep("^R[0-9]", names(data), value = TRUE)
    data$avg_neglogp <- rowMeans(-log10(data[, pval_cols] + 1e-10), na.rm = TRUE)
    
    # Get top N RBPs by significance
    top_rbps <- data %>%
      arrange(desc(avg_neglogp)) %>%
      head(input$rbp_top_n)
    
    # Define region labels based on event type
    if (input$rbp_event_type == "ex") {
      # 8 regions for exons
      region_labels <- c(
        "R1" = "Upstream\nExon\n(last 50bp)",
        "R2" = "Upstream\nIntron\n(first 50bp)",
        "R3" = "Upstream\nIntron\n(last 50bp)",
        "R4" = "Target\nExon\n(first 50bp)",
        "R5" = "Target\nExon\n(last 50bp)",
        "R6" = "Downstream\nIntron\n(first 50bp)",
        "R7" = "Downstream\nIntron\n(last 50bp)",
        "R8" = "Downstream\nExon\n(first 50bp)"
      )
    } else {
      # 5 regions for introns
      region_labels <- c(
        "R1" = "Upstream\nExon\n(first 50bp)",
        "R2" = "5' Splice Site\n(last 50bp\nof exon)",
        "R3" = "Intron\n(250bp)",
        "R4" = "3' Splice Site\n(first 50bp\nof exon)",
        "R5" = "Downstream\nExon\n(last 50bp)"
      )
    }
    
    # Create heatmap-style plot
    plot_data <- top_rbps %>%
      select(RBP_clean, all_of(pval_cols)) %>%
      pivot_longer(cols = all_of(pval_cols), names_to = "Region", values_to = "pvalue") %>%
      mutate(
        neglogp = -log10(pvalue + 1e-10),
        neglogp = pmin(neglogp, 4)  # Cap at 4
      )
    
    # Add region labels
    plot_data$Region_label <- region_labels[plot_data$Region]
    
    # Create ordered factors
    region_order <- names(region_labels)
    region_label_order <- unique(region_labels[region_order])
    
    # Order RBPs by average significance
    rbp_order <- top_rbps %>%
      arrange(desc(avg_neglogp)) %>%
      pull(RBP_clean) %>%
      unique()
    
    plot_data$RBP_clean <- factor(plot_data$RBP_clean, levels = rbp_order)
    plot_data$Region_label <- factor(plot_data$Region_label, levels = region_label_order)
    
    ggplot(plot_data, aes(x = Region_label, y = RBP_clean, fill = neglogp)) +
      geom_tile(color = "white", linewidth = 0.5) +
      scale_fill_gradient(low = "#E8E3F5", high = "#706993", name = "-log10(p)") +
      labs(
        title = paste("RBP Motif Enrichment -", input$rbp_drug, 
                      ifelse(input$rbp_event_type == "ex", "(Exons)", "(Introns)")),
        subtitle = paste("Top", input$rbp_top_n, "RBPs ranked by significance"),
        x = NULL,
        y = "RBP"
      ) +
      theme_publication() +
      theme(
        axis.text.x = element_text(angle = 0, size = 8, hjust = 0.5),
        axis.text.y = element_text(size = 9),
        legend.position = "right"
      )
  })
  
  output$download_rbp <- downloadHandler(
    filename = function() {
      paste0("rbp_analysis_", Sys.Date(), ".csv")
    },
    content = function(file) {
      data <- rbp_data()
      if (!is.null(data)) {
        write.csv(data, file, row.names = FALSE)
      }
    }
  )
  
  # ============================================================================
  # DATA TABLES TAB
  # ============================================================================
  
  table_filtered_data <- reactive({
    req(input$table_drug)
    
    data <- filtered_data()[[input$table_drug]]
    req(data)
    
    if (input$table_filter == "significant") {
      data <- filter(data, significant == TRUE)
    } else if (input$table_filter == "nonsig") {
      data <- filter(data, significant == FALSE)
    }
    
    # Select relevant columns for display
    display_cols <- c("GENE", "EVENT", "COORD", "LENGTH", "event_type", 
                      "deltapsi", "FDR", "direction")
    
    data <- data %>%
      select(any_of(display_cols))
    
    # Rename columns explicitly
    col_mapping <- c(
      "GENE" = "Gene",
      "EVENT" = "Event ID",
      "COORD" = "Coordinates",
      "LENGTH" = "Length (bp)",
      "event_type" = "Event Type",
      "deltapsi" = "ΔPSI",
      "FDR" = "FDR",
      "direction" = "Direction"
    )
    
    # Apply renaming for columns that exist
    existing_cols <- intersect(names(data), names(col_mapping))
    for (old_name in existing_cols) {
      names(data)[names(data) == old_name] <- col_mapping[old_name]
    }
    
    data
  })
  
  output$data_table <- renderDT({
    data <- table_filtered_data()
    
    datatable(
      data,
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel')
      ),
      rownames = FALSE,
      filter = 'top'
    ) %>%
      formatRound(columns = c("ΔPSI"), digits = 3) %>%
      formatSignif(columns = c("FDR"), digits = 3) %>%
      formatStyle(
        'Direction',
        backgroundColor = styleEqual(
          c('Increased', 'Decreased'),
          c(paste0(SPLICING_COLORS["Included"], "40"),
            paste0(SPLICING_COLORS["Skipped"], "40"))
        )
      )
  })
  
  # ============================================================================
  # DOWNLOAD HANDLERS
  # ============================================================================
  
  # Download volcano plot
  output$download_volcano <- downloadHandler(
    filename = function() {
      paste0("volcano_", tolower(gsub(" ", "_", input$volcano_drug)), "_", 
             Sys.Date(), ".pdf")
    },
    content = function(file) {
      data <- volcano_plot_data()
      
      color_map <- c(
        "Increased" = SPLICING_COLORS["Included"],
        "Decreased" = SPLICING_COLORS["Skipped"],
        "Not Significant" = SPLICING_COLORS["Unchanged"]
      )
      
      p <- ggplot(data, aes(x = deltapsi, y = negLog10FDR, color = color_group)) +
        geom_point(aes(size = negLog10FDR), alpha = 0.6) +
        scale_color_manual(values = color_map, name = "") +
        geom_hline(yintercept = -log10(input$fdr_threshold), 
                   linetype = "dashed", color = "gray50") +
        geom_vline(xintercept = c(-input$deltapsi_threshold, input$deltapsi_threshold), 
                   linetype = "dashed", color = "gray50") +
        labs(
          title = paste("Volcano Plot -", input$volcano_drug),
          x = "ΔPSI",
          y = "-log10(FDR)"
        ) +
        theme_publication()
      
      if (input$volcano_labels) {
        label_data <- filter(data, show_label)
        p <- p + geom_text_repel(
          data = label_data,
          aes(label = GENE),
          size = 3,
          max.overlaps = 20
        )
      }
      
      ggsave(file, p, width = 10, height = 8, dpi = 300)
    }
  )
  
  # Download table as CSV
  output$download_table_csv <- downloadHandler(
    filename = function() {
      paste0("splicing_events_", tolower(gsub(" ", "_", input$table_drug)), "_",
             Sys.Date(), ".csv")
    },
    content = function(file) {
      write_csv(table_filtered_data(), file)
    }
  )
  
  # Download table as Excel
  output$download_table_excel <- downloadHandler(
    filename = function() {
      paste0("splicing_events_", tolower(gsub(" ", "_", input$table_drug)), "_",
             Sys.Date(), ".xlsx")
    },
    content = function(file) {
      writexl::write_xlsx(table_filtered_data(), file)
    }
  )
  
  # Download all FDR tables
  output$download_all_fdr <- downloadHandler(
    filename = function() {
      paste0("all_fdr_tables_", Sys.Date(), ".zip")
    },
    content = function(file) {
      # Create temp directory
      temp_dir <- tempdir()
      
      # Save all tables
      for (drug in names(datasets)) {
        if (!is.null(datasets[[drug]])) {
          write_csv(datasets[[drug]], file.path(temp_dir, paste0(drug, "_fdr.csv")))
        }
      }
      
      # Create zip
      zip(file, list.files(temp_dir, pattern = "*_fdr.csv$", full.names = TRUE))
    }
  )
  
  # Download events barplot
  output$download_events <- downloadHandler(
    filename = function() { 
      paste0("events_barplot_", Sys.Date(), ".pdf") 
    },
    content = function(file) { 
      req(significant_data())
      
      event_data <- lapply(names(significant_data()), function(drug) {
        data <- significant_data()[[drug]]
        data %>%
          mutate(Drug = drug) %>%
          count(Drug, event_type)
      }) %>% bind_rows()
      
      p <- ggplot(event_data, aes(x = event_type, y = n, fill = Drug)) +
        geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.3) +
        geom_text(aes(label = n), position = position_dodge(width = 0.8), vjust = -0.5, size = 3) +
        scale_fill_manual(values = DRUG_COLORS) +
        labs(
          title = "Event Type Distribution by Drug",
          x = "Event Type",
          y = "Count"
        ) +
        theme_publication() +
        scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
      
      ggsave(file, p, width = 10, height = 6, dpi = 300)
    }
  )
  
  # Download upset plot
  output$download_upset <- downloadHandler(
    filename = function() { 
      paste0("upset_plot_", input$upset_level, "_", Sys.Date(), ".pdf") 
    },
    content = function(file) {
      req(upset_data())
      
      # Save as PDF with proper dimensions
      pdf(file, width = 12, height = 8)
      
      sets <- upset_data()
      if (length(sets) >= 2) {
        m <- make_comb_mat(sets)
        
        # Set colors
        comb_colors <- sapply(comb_name(m), function(comb_name) {
          sets_in_comb <- extract_comb(m, comb_name)
          sets_present <- names(sets_in_comb)[sets_in_comb]
          
          if (length(sets_present) == 1) {
            return(DRUG_COLORS[sets_present])
          } else if (length(sets_present) == 2) {
            return("grey40")
          } else {
            return("grey20")
          }
        })
        
        upset_plot <- UpSet(
          m,
          set_order = names(sets),
          comb_order = order(-comb_size(m)),
          top_annotation = upset_top_annotation(m, add_numbers = TRUE),
          left_annotation = upset_left_annotation(m, add_numbers = TRUE),
          right_annotation = NULL,
          comb_col = comb_colors,
          bg_col = c("grey95", "white")
        )
        
        draw(upset_plot)
      }
      
      dev.off()
    }
  )
  
  # Download heatmap (save as PNG since it's plotly)
  output$download_heatmap <- downloadHandler(
    filename = function() { 
      paste0("heatmap_", tolower(gsub(" ", "_", input$heatmap_drug)), "_", Sys.Date(), ".pdf") 
    },
    content = function(file) {
      req(input$heatmap_drug, input$heatmap_event_type, significant_data())
      
      # Note: This is a placeholder since plotly heatmaps are interactive
      # User should screenshot or export from the browser
      # For PDF, we'd need to recreate with ggplot2 or ComplexHeatmap
      pdf(file, width = 12, height = 10)
      plot.new()
      text(0.5, 0.5, paste("Interactive heatmap for", input$heatmap_drug, 
                           "\nPlease use the browser export function for full quality"),
           cex = 1.2)
      dev.off()
    }
  )
  
  # Download GO enrichment plot
  output$download_go <- downloadHandler(
    filename = function() { 
      paste0("go_enrichment_", tolower(gsub(" ", "_", input$go_drug)), "_", 
             input$go_ontology, "_", Sys.Date(), ".pdf") 
    },
    content = function(file) {
      data <- go_filtered_data()
      req(data)
      
      top_data <- data %>% head(input$go_top_n)
      
      p <- ggplot(top_data, aes(x = reorder(Description, -log10(p.adjust)), 
                                 y = -log10(p.adjust))) +
        geom_col(fill = DRUG_COLORS[input$go_drug], color = "black", linewidth = 0.3) +
        coord_flip() +
        labs(
          title = paste("GO Enrichment -", input$go_drug, "-", input$go_ontology),
          x = NULL,
          y = "-log10(Adjusted P-value)"
        ) +
        theme_publication() +
        scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
      
      ggsave(file, p, width = 10, height = 8, dpi = 300)
    }
  )
  
  # Download GC content plot
  output$download_gc <- downloadHandler(
    filename = function() { 
      paste0("gc_content_", input$gc_event_type, "_", Sys.Date(), ".pdf") 
    },
    content = function(file) {
      data <- gc_plot_data()
      req(data, nrow(data) > 0)
      
      splicing_colors <- c(
        "Retained" = "#D4A574",
        "Removed" = "#B85450",
        "Unchanged" = "#CCCCCC"
      )
      
      comparisons <- list(
        c("Retained", "Removed"),
        c("Retained", "Unchanged"),
        c("Removed", "Unchanged")
      )
      
      p <- ggplot(data, aes(x = splicing, y = GC, fill = splicing)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6, linewidth = 0.5) +
        geom_jitter(width = 0.15, alpha = 0.5, size = 1.2) +
        stat_compare_means(
          method = "wilcox.test",
          comparisons = comparisons,
          label = "p.signif",
          label.size = 4
        ) +
        scale_fill_manual(values = splicing_colors) +
        facet_wrap(~ Drug, nrow = 1) +
        labs(
          title = paste("GC Content Distribution -", 
                        ifelse(input$gc_event_type == "exon", "Exons", "Introns")),
          x = "Splicing Outcome",
          y = "GC Content (%)"
        ) +
        theme_publication() +
        theme(
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = rel(0.85)),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = rel(0.95)),
          panel.spacing = unit(1, "lines")
        )
      
      ggsave(file, p, width = 12, height = 5, dpi = 300)
    }
  )
  
  # Download length distribution plot
  output$download_length <- downloadHandler(
    filename = function() { 
      paste0("length_distribution_", input$length_event_type, "_", Sys.Date(), ".pdf") 
    },
    content = function(file) {
      data <- gc_plot_data()
      req(data, nrow(data) > 0)
      
      splicing_colors <- c(
        "Retained" = "#D4A574",
        "Removed" = "#B85450",
        "Unchanged" = "#CCCCCC"
      )
      
      comparisons <- list(
        c("Retained", "Removed"),
        c("Retained", "Unchanged"),
        c("Removed", "Unchanged")
      )
      
      threshold_95 <- quantile(data$LE_o, 0.95, na.rm = TRUE)
      data <- data %>%
        mutate(
          LE_o_display = pmin(LE_o, threshold_95),
          is_outlier = LE_o > threshold_95
        )
      
      p <- ggplot(data, aes(x = splicing, y = LE_o_display, fill = splicing)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6, linewidth = 0.5) +
        geom_jitter(width = 0.15, alpha = 0.5, size = 1.2) +
        stat_compare_means(
          method = "wilcox.test",
          comparisons = comparisons,
          label = "p.signif",
          label.size = 4
        ) +
        scale_fill_manual(values = splicing_colors) +
        facet_wrap(~ Drug, nrow = 1) +
        labs(
          title = paste("Length Distribution -", 
                        ifelse(input$length_event_type == "exon", "Exons", "Introns")),
          x = "Splicing Outcome",
          y = "Length (nt)"
        ) +
        theme_publication() +
        theme(
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = rel(0.85)),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = rel(0.95)),
          panel.spacing = unit(1, "lines")
        )
      
      n_outliers <- sum(data$is_outlier, na.rm = TRUE)
      if (n_outliers > 0) {
        p <- p + labs(caption = paste0("Note: ", n_outliers, " outliers > ", 
                                        round(threshold_95), " nt clipped for display"))
      }
      
      ggsave(file, p, width = 12, height = 5, dpi = 300)
    }
  )
  
  # Download RBP analysis
  output$download_rbp <- downloadHandler(
    filename = function() { 
      paste0("rbp_analysis_", input$rbp_event_type, "_", Sys.Date(), ".pdf") 
    },
    content = function(file) {
      p <- rbp_plot()
      req(p)
      ggsave(file, p, width = 12, height = 8, dpi = 300)
    }
  )
  
  # Download all GO results as ZIP
  output$download_all_go <- downloadHandler(
    filename = function() { 
      paste0("all_go_results_", Sys.Date(), ".zip") 
    },
    content = function(file) {
      temp_dir <- tempdir()
      zip_files <- c()
      
      for (drug_name in names(go_data)) {
        for (ont in names(go_data[[drug_name]])) {
          go_result <- go_data[[drug_name]][[ont]]
          if (!is.null(go_result) && nrow(go_result) > 0) {
            csv_file <- file.path(temp_dir, paste0(drug_name, "_GO_", ont, ".csv"))
            write_csv(go_result, csv_file)
            zip_files <- c(zip_files, csv_file)
          }
        }
      }
      
      if (length(zip_files) > 0) {
        zip::zip(file, zip_files, mode = "cherry-pick")
      }
    }
  )
}
