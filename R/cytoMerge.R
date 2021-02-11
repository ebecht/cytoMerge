## library(parallel) ## pop when packaging

## library(flowCore)
## library(uwot)
## library(preprocessCore)
## library(xgboost)
## library(tools)
## options("bitmapType" = "cairo")

## ##################
## Arguments
## ##################

## train_set_fraction = 0.5

## input_files = list.files("/home/etienne/Desktop/test2/input/", full.names = TRUE, pattern = ".fcs")
## ## excluded_markers = c("Viability", "Time")
## excluded_markers = c("191Ir_DNA1", "193Ir_DNA2", "195Pt", "196Pt", "SampleID")
## rds_dir = tempdir()
## output_dir = "/home/etienne/Desktop/test2/out"
## writeTmpFilesToDisk = TRUE
## verbose = TRUE
## plotResults = TRUE

## xgboost_params = list(
##     nrounds = 100,
##     eta = 0.05,
##     verbose = 0L,
##     nthread = ifelse(interactive(), min(6L, detectCores()), detectCores())
## )
## umap_params = list(
##     n_neighbors = 50L,
##     n_epochs = 1000L,
##     verbose = FALSE,
##     n_threads = detectCores(),
##     n_sgd_threads = detectCores()
## )
## chop_quantiles = 0.005
## do_UMAP_and_plot = TRUE

#' Merge partially overlapping flow cytometry panels from multiple files that profile the same biological sample. This function first logicle-transform the data, then applies quantile-normalization across the samples. It then uses XGBoost to predict all markers across all files (including markers that are overlapping - these can be used as a prediction quality control).
#' @export
#' @param input_files Character vector of input FCS files: from the same biological samples, but with distinct and partially-overlapping antibody panels
#' @param excluded_markers Markers to ignore (according to the FCS files' "name" parameter description - or "desc" if "name" is missing). For instance, Viability, Time, DNA content ...
#' @param rds_dir Directory to store temporary (.Rds) files such as the regression models or normalized and harmonized data matrices. Defaults to a temporary subdirectory.
#' @param output_dir Directory to save the output. Must be empty at the start of execution. Defaults to a temporary subdirectory.
#' @param writeTmpFilesToDisk Boolean. Whether to save .Rds intermediary files to disk for further analyses. Defaults to FALSE.
#' @param verbose Boolean. Whether to plot progress information.
#' @param do_UMAP_and_plot Boolean. Whether to compute UMAP dimensionality reduction of the backbone and prediction-enriched expression matrices, and to overlay measured and predicted markers' expression on the UMAP embeddings. Setting this to TRUE will produce PDF files in the output directory (in addition to the FCS files). Setting this to TRUE will also include UMAP embeddings into the exported FCS files.
#' @param xgboost_params List of named arguments passed to xgboost::xgboost().
#' @param umap_params List of named arguments passed to uwot::umap().
#' @param chop_quantiles numeric (should be close to 0). If do_UMAP_and_plot is TRUE, on the PDF plots this will clip the top and bottom parts of the data to the corresponding quantile. Defaults to 0.005 (so the bottom 0.005 and top 0.995 quantiles are clipped). This only affect the color-mapping of the plots, not the predictions or UMAP embeddings.
#' @param train_set_fraction Fraction of the data to use as training set. Performance will be computed on the train and test sets. Defaults to 0.8.
#' @importFrom flowCore read.FCS
#' @importFrom flowCore estimateLogicle
#' @importFrom flowCore write.FCS
#' @importFrom flowCore inverseLogicleTransform
#' @importFrom flowCore fr_append_cols
#' @importFrom uwot umap
#' @importFrom xgboost xgboost
#' @importFrom tools file_path_sans_ext
#' @importFrom grDevices hcl.colors
#' @importFrom graphics abline axis image legend mtext title
#' @importFrom stats approx cor predict quantile runif setNames
#' @importFrom utils combn write.csv
#' @importFrom Biobase pData
#' @importFrom Biobase pData<-
#' @importFrom Biobase exprs
merge_flow_panels = function(
                             input_files,
                             excluded_markers = character(),
                             rds_dir = paste0(tempdir(), "tmp"),
                             output_dir = paste0(tempdir(), "out"),
                             writeTmpFilesToDisk = FALSE,
                             verbose = TRUE,
                             xgboost_params = list(
                                 nrounds = 100,
                                 eta = 0.05,
                                 verbose = 0L
                             ),
                             umap_params = list(
                                 n_neighbors = 50L,
                                 n_epochs = 1000L,
                                 verbose = FALSE
                             ),
                             chop_quantiles = 0.005,
                             do_UMAP_and_plot = TRUE,
                             train_set_fraction = 0.8
                             )
{

    ## ##################
    ## Sanity checks
    ## ##################
    w = duplicated(basename(input_files))
    if(any(w)){
        stop("input files cannot have the same filenames : ", basename(input_files)[input_files %in% input_files[w]])
    }
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    if(length(list.files(output_dir) > 0)){
        stop("Output directory ", output_dir, " is not empty")
    }

    ## ##################
    ## Data import and parsing
    ## ##################

    flowFrames_list = lapply(
        input_files,
        function(file){
            res = read.FCS(file, truncate_max_range = FALSE)
            pData(res@parameters)$desc[is.na(pData(res@parameters)$desc)] = pData(res@parameters)$name[is.na(pData(res@parameters)$desc)]
            res
        }
    )
    names(flowFrames_list) = input_files

    panels_forward = lapply(
        flowFrames_list,
        function(x){
            res = setNames(
                pData(x@parameters)$desc,
                pData(x@parameters)$name
            )
            res[!res %in% excluded_markers]
        }
    )
    panels_rev = lapply(panels_forward, function(x){setNames(names(x), x)})

    n_events = sapply(flowFrames_list, nrow, simplify = FALSE)

    files_pairs = combn(names(flowFrames_list), 2, simplify = FALSE)

    ## Normalize
    ## Define train set
    ## Train models
    ## Do predictions
    ## Do UMAPs
    ## Export
    ## Do plots

    ## ##################
    ## Define training set 
    ## ##################

    train_sets = lapply(
        flowFrames_list,
        function(x){
            w = rep(FALSE, nrow(x))
            n_train = ceiling(train_set_fraction * nrow(x))
            if(n_train > nrow(x)){
                stop("train_set_fraction is too large, should be (much) lower than 1")
            }
            w[sample(seq_len(nrow(x)), n_train)] = TRUE
            w
        }
    )
    if(writeTmpFilesToDisk){
        saveRDS(train_sets, file = file.path(rds_dir, "train_sets.Rds"))
    }
    
    ## ##################
    ## Predictions
    ## ##################
    lgcl_list = lapply(
        flowFrames_list,
        function(ff){
            estimateLogicle(ff, colnames(exprs(ff[, !ff@parameters$desc %in% excluded_markers])))
        }
    )
    if(writeTmpFilesToDisk){
        saveRDS(lgcl_list, file = file.path(rds_dir, "lgcl_list.Rds"))
    }

    predictions = lapply(
        files_pairs,
        function(files_pair){## Directories
            rds_subdir = file.path(rds_dir, paste0(tools::file_path_sans_ext(basename(files_pair)), collapse = " / "))
            dir.create(rds_subdir, showWarnings = FALSE, recursive = TRUE)

            ## Definition of predictor variables
            predictors = panels_rev[files_pair]
            predictors = lapply(
                predictors,
                "[",
                setdiff(sort(Reduce(intersect, lapply(predictors, names))), excluded_markers)
            )
            if(writeTmpFilesToDisk){
                saveRDS(predictors, file = file.path(rds_subdir, "predictors.Rds"))
            }
            
            ## Definition of target variables
            targets = panels_rev[files_pair]
            targets = sapply(
                names(targets),
                function(x){
                    panels_rev[[x]][
                        panels_forward[[x]][setdiff(
                                          targets[[x]],
                                          predictors[[x]]                            
                                      )]
                    ]
                },
                simplify = FALSE
            )
            if(writeTmpFilesToDisk){
                saveRDS(targets, file = file.path(rds_subdir, "targets.Rds"))
            }
            
            fs = lapply(
                files_pair,
                function(x){
                    ff = flowFrames_list[[x]]
                    ff[, c(predictors[[x]], targets[[x]])]
                }
            )
            names(fs) = files_pair
            
            xp = sapply(
                names(fs),
                function(id){
                    exprs(flowCore::transform(fs[[id]], lgcl_list[[id]]))
                },
                simplify = FALSE
            )
            if(writeTmpFilesToDisk){
                saveRDS(xp, file.path(rds_subdir, "xp_transformed.Rds"))
            }
            
            xp_transformed_backbone = sapply(
                names(xp),
                function(id){
                    xp[[id]][, predictors[[id]]]
                },
                simplify = FALSE
            )
            saveRDS(xp, file.path(rds_subdir, "xp_transformed_backbone.Rds"))

            xp_transformed_targets = sapply(
                names(xp),
                function(id){
                    xp[[id]][, targets[[id]]]
                },
                simplify = FALSE
            )
            if(writeTmpFilesToDisk){
                saveRDS(xp, file.path(rds_subdir, "xp_transformed_targets.Rds"))
            }
            
            ## Normalization of the backbone: Z-score + quantile normalization for each column across FCS files
            xp_normalized_backbone = lapply(xp_transformed_backbone, scale)
            events.code = rep(names(xp_normalized_backbone), times = sapply(xp_normalized_backbone, nrow))
            xp_normalized_backbone = do.call(rbind, xp_normalized_backbone)
            xp_normalized_backbone = QN.xp(xp_normalized_backbone, events.code)       
            xp_normalized_backbone = lapply(split(as.data.frame(xp_normalized_backbone), events.code), as.matrix)
            for(id in names(xp_normalized_backbone)){
                colnames(xp_normalized_backbone[[id]]) = predictors[[id]]
            }
            if(writeTmpFilesToDisk){
                saveRDS(xp_normalized_backbone, file.path(rds_subdir, "xp_normalized_backbone.Rds"))
                saveRDS(events.code, file.path(rds_subdir, "events.code.Rds"))
            }

            ## Training models for leave-one-out backbone prediction
            if(verbose){
                cat("\n", "Training leave-one-out backbone models", "\n")
            }
            models_leave_one_out_backbone = sapply(
                names(targets),
                function(id){
                    if(verbose){
                        cat("\n\t", tools::file_path_sans_ext(basename(id)), "\n")
                    }
                    sapply(
                        unname(predictors[[id]]),
                        function(predictor){
                            if(verbose){
                                cat("\t\t", predictor, " (", panels_forward[[id]][predictor], ")", "\n", sep = "")
                            }
                            do.call(
                                xgboost,
                                c(
                                    list(
                                        data = xp_normalized_backbone[[id]][train_sets[[id]], setdiff(colnames(xp_normalized_backbone[[id]]), predictor)],
                                        label = xp_normalized_backbone[[id]][train_sets[[id]], predictor]
                                    ),
                                    xgboost_params
                                )
                            )
                        },
                        simplify = FALSE
                    )
                },
                simplify = FALSE
            )
            if(writeTmpFilesToDisk){
                saveRDS(models_leave_one_out_backbone, file = file.path(rds_subdir, "models_leave_one_out_backbone.Rds"))
            }
            
            ## Training models for targets (NB: for core we are actually training the models multiple times... not computationnally efficient and redudant)
            if(verbose){
                cat("\n", "Training variable-targets models", "\n")
            }
            models_variable_targets = sapply(
                names(targets),
                function(id){
                    if(verbose){
                        cat("\n\t", tools::file_path_sans_ext(basename(id)), "\n")
                    }
                    sapply(
                        unname(targets[[id]]),
                        function(target){
                            if(verbose){
                                cat("\t\t", target, " (", panels_forward[[id]][target], ")", "\n", sep = "")
                            }
                            do.call(
                                xgboost,
                                c(
                                    list(
                                        data = xp_normalized_backbone[[id]][train_sets[[id]], ],
                                        label = xp_transformed_targets[[id]][train_sets[[id]], target]
                                    ),
                                    xgboost_params
                                )
                            )
                        },
                        simplify = FALSE
                    )
                },
                simplify = FALSE
            )
            if(writeTmpFilesToDisk){
                saveRDS(models_variable_targets, file = file.path(rds_subdir, "models_variable_targets.Rds"))
            }

            if(verbose){
                cat("\n", "Predicting from leave-one-out models", "\n")
            }

            ## Prediction for core panel
            predictions_leave_one_out_backbone = sapply(
                names(models_leave_one_out_backbone),
                function(id){
                    if(verbose){
                        cat("\n\t", tools::file_path_sans_ext(basename(id)), "\n")
                    }
                    sapply(
                        names(models_leave_one_out_backbone[[id]]),
                        function(channel){
                            if(verbose){
                                cat("\t\t", channel," (", panels_forward[[id]][channel], ")", "\n", sep = "")
                            }
                            cn = colnames(xp_normalized_backbone[[id]])
                            predictors_matrix = xp_normalized_backbone
                            predictors_matrix = lapply(predictors_matrix, get("colnames<-"), value = cn)
                            predictors_matrix = do.call(rbind, predictors_matrix)
                            measured = predictors_matrix[, channel]
                            predictors_matrix = predictors_matrix[, setdiff(cn, channel)]
                            model = models_leave_one_out_backbone[[id]][[channel]]
                            predictions = predict(model, predictors_matrix)
                            train_set = rep(0, length(predictions))
                            train_set[events.code == id] = train_sets[[id]]
                            cbind(
                                measured = measured,
                                predicted = predictions,
                                train_set = train_set
                            )
                        },
                        simplify = FALSE
                    )
                },
                simplify = FALSE
            )
            if(writeTmpFilesToDisk){
                saveRDS(predictions_leave_one_out_backbone, file = file.path(rds_subdir, "predictions_leave_one_out_backbone.Rds"))
            }
            
            if(verbose){
                cat("\n", "Predicting from variable-targets models", "\n")
            }
            predictions_variable_targets = sapply(
                names(models_variable_targets),
                function(id){
                    if(verbose){
                        cat("\n\t", tools::file_path_sans_ext(basename(id)), "\n")
                    }
                    sapply(
                        names(models_variable_targets[[id]]),
                        function(channel){
                            if(verbose){
                                cat("\t\t", channel," (", panels_forward[[id]][channel], ")", "\n", sep = "")
                            }
                            cn = colnames(xp_normalized_backbone[[id]])
                            predictors_matrix = xp_normalized_backbone
                            predictors_matrix = lapply(predictors_matrix, get("colnames<-"), value = cn)
                            predictors_matrix = do.call(rbind, predictors_matrix)
                            measured = xp_transformed_targets[[id]][, channel]
                            model = models_variable_targets[[id]][[channel]]
                            predictions = predict(model, predictors_matrix[, model$feature_names])
                            train_set = rep(0, length(predictions))
                            train_set[events.code == id] = train_sets[[id]]
                            result = cbind(
                                predicted = predictions,
                                train_set = train_set
                            )
                            result = cbind(
                                measured = NA,
                                result
                            )
                            result[events.code == id, "measured"] = measured
                            result
                        },
                        simplify = FALSE
                    )
                },
                simplify = FALSE
            )
            if(writeTmpFilesToDisk){
                saveRDS(predictions_variable_targets, file = file.path(rds_subdir, "predictions_variable_targets.Rds"))
            }
            cat("\n")
            return(
                list(
                    predictions_leave_one_out_backbone = predictions_leave_one_out_backbone,
                    predictions_variable_targets = predictions_variable_targets,
                    targets = targets,
                    predictors = predictors,
                    events.code = events.code
                )
            )
        }
    )

    ## Tie everything together: output FCS files that for each event outputs all the possible predictions (for all [panels x instruments] output all predicted markers
    ## One prediction = panel x instrument x target x panels_pair_including_the_panel (a "backbone")
    ## So the number of predicted targets varies because sometimes it is actually measured on a subset of the dataset. In this case we just use the prediction from the leave one out model?
    ## In this case, predictions are each marker from each panel predicted
    if(verbose){
        cat("\n", "Consolidating predictions", "\n")
    }
    events.code_sub = rep(input_files, unlist(n_events))

    all_columns = do.call(
        c,
        unname(
            sapply(
                names(panels_rev),
                function(x){
                    paste0(unname(panels_rev[[x]]), "_predicted_from_", x)
                },
                simplify = FALSE
            )
        )
    )

    output_matrix = matrix(
        ncol = length(all_columns),
        nrow = length(events.code_sub),
        data = NA,
        dimnames = list(
            NULL,
            all_columns
        )
    )

    for(i in seq_len(length(files_pairs))){
        files_pair = files_pairs[[i]]
        
        ## Verbose progress
        if(verbose){
            cat("\n\t", paste0(tools::file_path_sans_ext(basename(files_pair)), collapse = " / "))
        }
        
        ## Directories
        rds_subdir = file.path(rds_dir, paste0(tools::file_path_sans_ext(basename(files_pair)), collapse = " / "))

        ## Load preprocessed data
        current_environment = environment()
        invisible(
            sapply(
                names(predictions[[i]]),
                function(object_name){
                    assign(
                        object_name,
                        predictions[[i]][[object_name]],
                        envir = current_environment
                    )
                }
            )
        )

        ## Concatenate predictions for leave-one-out backbone
        predictions_leave_one_out_backbone_concatenated = do.call(
            cbind,
            sapply(
                names(predictions_leave_one_out_backbone),
                function(x){
                    res = do.call(
                        cbind,
                        sapply(
                            predictions_leave_one_out_backbone[[x]],
                            function(y){
                                y[, "predicted"]
                            },
                            simplify = FALSE
                        )
                    )
                    colnames(res) = paste0(colnames(res), "_predicted_from_", x)
                    res
                },
                simplify = FALSE
            )
        )
        if(writeTmpFilesToDisk){
            saveRDS(predictions_leave_one_out_backbone, file = file.path(rds_subdir, "predictions_leave_one_out_backbone.Rds"))
        }
        
        ## Concatenate predictions for variable targets
        predictions_variable_targets_concatenated = do.call(
            cbind,
            sapply(
                names(predictions_variable_targets),
                function(x){
                    res = do.call(
                        cbind,
                        sapply(
                            predictions_variable_targets[[x]],
                            function(y){
                                y[, "predicted"]
                            },
                            simplify = FALSE
                        )
                    )
                    colnames(res) = paste0(colnames(res), "_predicted_from_", x)
                    res
                },
                simplify = FALSE
            )
        )
        if(writeTmpFilesToDisk){
            saveRDS(predictions_variable_targets_concatenated, file = file.path(rds_subdir, "predictions_variable_targets_concatenated.Rds"))
        }
        
        for(id in unique(events.code)){
            w_total = events.code_sub == id
            w_pred = events.code == id
            output_matrix[w_total, colnames(predictions_variable_targets_concatenated)] = predictions_variable_targets_concatenated[w_pred, ]
            output_matrix[w_total, colnames(predictions_leave_one_out_backbone_concatenated)] = predictions_leave_one_out_backbone_concatenated[w_pred, ]
        }
    }
    if(writeTmpFilesToDisk){
        saveRDS(output_matrix, file = file.path(rds_dir, "output_matrix.Rds"))
        saveRDS(events.code_sub, file = file.path(rds_dir, "events.code_sub.Rds"))
    }
    if(verbose){
        cat("\n")
    }

    if(verbose){
        cat("\n", "Exporting FCS files and plotting", "\n")
    }
    output_matrix_raw = output_matrix
    for(cn in colnames(output_matrix)){
        id = strsplit(cn, "_predicted_from_")[[1]]
        channel = id[[1]]
        id = id[[2]]
        ilgcl = inverseLogicleTransform(lgcl_list[[id]])
        output_matrix[, cn] = ilgcl@transforms[[channel]]@f(output_matrix[, cn])
    }
    output_matrix = split_matrix(output_matrix, events.code_sub)
    output_matrix_raw = split_matrix(output_matrix_raw, events.code_sub)

    fcs_files = lapply(
        names(output_matrix),
        function(id){
            if(verbose){
                cat("\n\t", tools::file_path_sans_ext(basename(id)), "\n")
            }
            
            input_file = id
            fcs = read.FCS(
                id,
                truncate_max_range = FALSE
            )
            nchannels = ncol(fcs)
            fcs_out = fr_append_cols(fcs, output_matrix[[id]])
            names = fcs_out@parameters$name[(nchannels+1):ncol(fcs_out)]
            files = strsplit(names, "_predicted_from_")
            channels = sapply(files, "[", 1)
            files = sapply(files, "[", 2)
            desc = sapply(
                seq_along(names),
                function(i){
                    paste0(panels_forward[[files[i]]][channels[i]], "_predicted_from_", tools::file_path_sans_ext(basename(files[i])))
                }
            )
            fcs_out@parameters$desc = c(fcs_out@parameters$desc[1:nchannels], desc)
            fcs_out@parameters$desc[is.na(fcs_out@parameters$desc)] = fcs_out@parameters$name[is.na(fcs_out@parameters$desc)]

            if(do_UMAP_and_plot){
                chans_sub = panels_forward[[id]][names(lgcl_list[[id]]@transforms)]
                chans_sub = chans_sub[chans_sub != "Viability"]
                umap_backbone = fcs_out
                umap_backbone = exprs(flowCore::transform(umap_backbone[, 1:nchannels], lgcl_list[[id]]))
                umap_backbone = umap_backbone[, names(chans_sub)]
                
                umap_backbone = do.call(uwot::umap, c(list(X = umap_backbone), umap_params))
                colnames(umap_backbone) = paste0("UMAP", 1:2, "_backbone")
                
                umap_augmented = fcs_out
                umap_augmented = exprs(flowCore::transform(umap_augmented[, 1:nchannels], lgcl_list[[id]]))
                umap_augmented = umap_augmented[, names(chans_sub)]
                umap_augmented = cbind(umap_augmented, output_matrix_raw[[id]])
                umap_augmented = do.call(uwot::umap, c(list(X = umap_augmented, init = umap_backbone), umap_params))
                colnames(umap_augmented) = paste0("UMAP", 1:2, "_augmented")

                fcs_out = fr_append_cols(fcs_out, cbind(umap_backbone, umap_augmented, training_set = as.numeric(train_sets[[id]])))
                
                all_measured = exprs(flowCore::transform(fcs_out[, 1:nchannels], lgcl_list[[id]]))[, names(chans_sub)]
                colnames(all_measured) = paste0(panels_forward[[id]][colnames(all_measured)], "_measured")
                all_pred = output_matrix_raw[[id]]
                colnames(all_pred) = setNames(fcs_out@parameters$desc, fcs_out@parameters$name)[colnames(all_pred)]
                all_pred = apply(
                    all_pred,
                    2,
                    function(x){
                        q = quantile(x, c(chop_quantiles, 1 - chop_quantiles))
                        x[x<q[1]] = q[1]
                        x[x>q[2]] = q[2]
                        x
                    }
                )
                
                all_data = cbind(
                    all_measured,
                    all_pred
                )
                all_data = all_data[, sort(colnames(all_data))]
                all_data = cbind(
                    all_data,
                    umap_backbone,
                    umap_augmented
                )
                
                for(type in c("backbone", "augmented")){
                    color_biplot_by_channels(
                                       all_data,
                                       x_axis = paste0("UMAP1_", type),
                                       y_axis = paste0("UMAP2_", type),
                                       pch = 16,
                                       cex = min(1, 1.1 - 0.15 * log10(nrow(all_data))),
                                       palette = hcl.colors(100, "viridis"),
                                       file_name = file.path(output_dir, paste0(sub(paste0(".", tools::file_ext(id)), "", basename(id), fixed = TRUE), "_umap_", type, "_colored_by_expression.pdf")),
                                       global_across_channels = FALSE,
                                       resolution = 150,
                                       raster.height = 1000,
                                       raster.width = 1000
                                   )
                }        
            }
            write.FCS(fcs_out, filename = file.path(output_dir, basename(id)))
            return(fcs_out)
        }
    )

    ## ##################
    ## Performance assesment
    ## ##################
    if(verbose){
        cat("\n", "Assessing performance and plotting", "\n")
    }
    for(i in seq_len(length(files_pairs))){
        files_pair = files_pairs[[i]]

        ## Verbose progress
        if(verbose){
            cat("\n\t", paste0(tools::file_path_sans_ext(basename(files_pair)), collapse = " / "), "\n")
        }
        
        ## Directories
        rds_subdir = file.path(rds_dir, paste0(tools::file_path_sans_ext(basename(files_pair)), collapse = " / "))
        
        ## Load preprocessed data
        current_environment = environment()
        invisible(
            sapply(
                names(predictions[[i]]),
                function(object_name){
                    assign(
                        object_name,
                        predictions[[i]][[object_name]],
                        envir = current_environment
                    )
                }
            )
        )

        all_predictions = sapply(
            names(predictions_leave_one_out_backbone),
            function(x){
                res = c(predictions_leave_one_out_backbone[[x]], predictions_variable_targets[[x]])
                res[sort(names(res))]
            },
            simplify = FALSE
        )

        perf = sapply(
            names(all_predictions),
            function(id){
                graph_dir = file.path(output_dir, "predictions_quality", paste0(tools::file_path_sans_ext(basename(files_pair)), collapse = " _ "), tools::file_path_sans_ext(basename(id)))
                cat("\t\t\t", tools::file_path_sans_ext(basename(id)), "\n")
                dir.create(graph_dir, showWarnings = FALSE, recursive = TRUE)
                sapply(
                    names(all_predictions[[id]]),
                    function(channel){
                        data = all_predictions[[id]][[channel]]
                        data = subset(data, !is.na(data[, "measured"]))
                        data = split_matrix(data[, c("measured", "predicted")], data[, "train_set"])
                        png(file.path(graph_dir, make.names(paste0(channel, "_", panels_forward[[id]][channel],".png"))), width = 4000, height = 2000, res = 300)
                        par(mfrow = c(1, 2), oma = c(2, 2, 2, 2))
                        x = data[["1"]][, "measured"]
                        y = data[["1"]][, "predicted"]
                        freqplot(x, y, main = "train set")
                        title(xlab = "measured", ylab = "predicted")
                        R2_train = cor(x, y)
                        legend(x = "topleft", legend = paste0("R2 = ", signif(R2_train^2, 2)), bty = "n")
                        abline(a = 0, b = 1, lty = 3)
                        main = "train set"
                        x = data[["0"]][, "measured"]
                        y = data[["0"]][, "predicted"]
                        freqplot(x, y, main = "test set")
                        title(xlab = "measured", ylab = "predicted")
                        R2_pred = cor(x, y)
                        main = "test set"
                        legend(x = "topleft", legend = paste0("R2 = ", signif(R2_pred^2, 2)), bty = "n")
                        abline(a = 0, b = 1, lty = 3)                            
                        mtext(side = 3, at = 0.5, text = paste0(channel, " (", panels_forward[[id]][channel], ")"), outer = TRUE)
                        dev.off()

                        return(c(R2_train = R2_train, R2_pred = R2_pred))
                    }
                )
            },
            simplify = FALSE
        )
        perf = lapply(perf, t)
        perf = lapply(perf, as.data.frame)
        perf = sapply(
            names(perf),
            function(x){
                cbind(
                    dataset = x,
                    channel = rownames(perf[[x]]),
                    target = panels_forward[[x]][rownames(perf[[x]])],
                    files_pair = paste0(tools::file_path_sans_ext(basename(files_pair)), collapse = " / "),
                    model_type = ifelse(rownames(perf[[x]]) %in% predictors[[x]], "leave_one_out_backbone", "cross_panel_targets"),
                    perf[[x]]
                )
            },
            simplify = FALSE
        )
        perf = do.call(rbind, perf)
        rownames(perf) = NULL
        write.csv(perf, file.path(output_dir, "predictions_quality", paste0(tools::file_path_sans_ext(basename(files_pair)), collapse = " _ "), "performance_summary_table.csv"), row.names = FALSE)

        NULL
    }

    cat("\nResults are available in ", output_dir, "\n")
}

## ##################
## Helper functions
## ##################

#' Quantile normalization for multiple vectors of possibly different sizes
#' @importFrom preprocessCore normalize.quantiles
#' @param lov a list of numeric vectors
#' @note Uses linear interpolation if the sizes are different
QN = function (lov){
    w.n = which.max(sapply(lov, length))
    n = length(lov[[w.n]])
    o = lapply(lov, order)
    or = lapply(o, order)
    lov = lapply(1:length(lov), function(x) lov[[x]][o[[x]]])
    lov.interpolated = lapply(lov, approx, n = n)
    lov.interpolated = lapply(lov.interpolated, function(x) x$y)
    lov.interpolated = do.call(cbind, lov.interpolated)
    res = normalize.quantiles(lov.interpolated)
    lapply(1:length(lov), function(i) {
        unname(quantile(res[, i], seq(0, 1, length.out = length(lov[[i]])))[or[[i]]])
    })
}

#' Quantile normalization on a matrix + vector of groups (splitting by rows) 
#' @param xp_mat a numeric matrix
#' @param groups a vector of length nrow(xp_mat)
QN.xp = function (xp_mat, groups){
    o = order(unlist(split(1:nrow(xp_mat), groups), use.names = FALSE))
    apply(xp_mat, 2, function(x) {
        unlist(QN(split(x, groups)))[o]
    })
}

#' splitting a matrix to a list by a margin according to a vector of groups
#' @param mat a numeric matrix
#' @param vector a vector of groups labels
#' @param byrow boolean. Whether to split by rows or columns
split_matrix = function (mat, vector, byrow = TRUE){
    if (byrow & nrow(mat) != length(vector)) {
        stop("if byrow=TRUE, vector's length should have length nrow(mat)")
    }
    else if (!byrow & ncol(mat) != length(vector)) {
        !byrow & ncol(mat) != length(vector)
        stop("if byrow=FALSE, vector's length should have length ncol(mat)")
    }
    if (byrow) {
        levels = split(1:nrow(mat), vector)
        res = lapply(levels, function(x) {
            mat[x, , drop = FALSE]
        })
    }
    else {
        levels = split(1:ncol(mat), vector)
        res = lapply(levels, function(x) {
            mat[, x, drop = FALSE]
        })
    }
    res
}

#' Heatmap for 2D histograms
#' @param x numerical vector
#' @param y numerical vector
#' @param breaks Number of bins
#' @param na.rm Whether to exclude NAs from x and y
#' @param palette color palette for the heatmap
#' @param add_white Boolean. Whether to add white at the beginning of the palette (for counts of 0 in the 2D histogram)
#' @param ... passed to image()
freqplot = function (x, y, breaks = 200, na.rm = TRUE, palette = rev(c("#A50026", 
    "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#FFFFBF", "#E0F3F8", 
    "#ABD9E9", "#74ADD1", "#4575B4", "#313695")), add_white = TRUE, 
    ...){
    w = is.na(x) | is.na(y) | is.nan(x) | is.nan(y) | is.infinite(x) | 
        is.infinite(y)
    if (any(w)) {
        if (na.rm) {
            x = x[!w]
            y = y[!w]
        }
        else {
            stop("NA values found and na.rm is FALSE")
        }
    }
    w.x = length(unique(x)) > 1
    w.y = length(unique(y)) > 1
    if (w.x) {
        breaks.x = seq(min(x), max(x), length.out = breaks)
        labels.x = breaks.x[-length(breaks.x)]
        X = cut(x, breaks = breaks.x, labels = labels.x, include.lowest = TRUE)
    }
    else {
        X = x
    }
    if (w.y) {
        breaks.y = seq(min(y), max(y), length.out = breaks)
        labels.y = breaks.y[-length(breaks.y)]
        Y = cut(y, breaks = breaks.y, labels = labels.y, include.lowest = TRUE)
    }
    else {
        Y = y
    }
    tab = log10(1 + table(X, Y))
    if (length(x) < 1 | length(y) < 1) {
        plot.new()
        return(tab)
    }
    if (w.x & w.y) {
        if (add_white) {
            null_color = "white"
        }
        else {
            null_color = NULL
        }
        image(tab, col = c(null_color, colorRampPalette(palette)(100)), 
            x = breaks.x, y = breaks.y, xaxt = "n", yaxt = "n", 
            xlab = "", ylab = "", bty = "n", ...)
        ticks = seq(0, 1, by = 0.25)
        if (par("xaxt") != "n") {
            axis(side = 1, at = quantile(breaks.x, ticks), labels = signif(quantile(breaks.x, 
                ticks), 2), line = 0.5)
        }
        if (par("yaxt") != "n") {
            axis(side = 2, at = quantile(breaks.y, ticks), labels = signif(quantile(breaks.y, 
                ticks), 2), line = 0.5)
        }
    }
    else {
        if (!w.x) {
            X = runif(length(x))
        }
        else {
            X = x
        }
        if (!w.y) {
            Y = runif(length(y))
        }
        else {
            Y = y
        }
        freqplot(X, Y, breaks = breaks, na.rm = na.rm, ...)
    }
    tab
}


#' Colors points of a biplot (2d-tSNE, 2d-PCA...) by the intensity of channels for each flowframe in the flowset
#' @param matrix A matrix
#' @param x_axis A column name of matrix used in biplots as the x axis
#' @param y_axis A column name of matrix used in biplots as the y axis
#' @param global_across_channels Boolean specificying whether the color key should be calibrated across all channels or individually per channel.
#' @param palette A vector of colors that'll be passed to colorRampPalette
#' @param resolution The resolution of the files that'll be imported in the pdf. Default is 72, increase for higher resolution. You may need to enlarge rasters accordingly.
#' @param data_transformation_reverse The colors will be linearly scaled across the range of the data. If the data was transformed, you may however want the labels which appear in your color-keys to reflect the raw intensities. In this case, this should be the inverse of the function you used to transform your data
#' @param file_name String of the location of the output file (should end in .pdf)
#' @param raster.width Width of each embedded raster. Defaults to 480
#' @param raster.height height of each embedded raster. Defaults to 480
#' @param ... passed to plot (suggested: pch=16, cex=0.5 or less)
#' @return NULL
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices png
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @importFrom raster tmpDir
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics plot.new
#' @importFrom graphics rect
#' @importFrom graphics text
#' @importFrom grid grid.raster
#' @importFrom grid grid.newpage
#' @importFrom grid grid.text
#' @importFrom grid unit
#' @importFrom grid gpar
#' @importFrom png readPNG
#' @importFrom utils tail
#' @note Since pdf files are vectorized, they can get really big if a lot of data point are plotted. This function thus used bitmap images that are stored in a temporary directory (tmpDir()) and then import them in a single pdf. If you're interested in using the bitmap images, you can fetch them in tmpDir()
#' @noRd

color_biplot_by_channels <- function(
                                     matrix,
                                     x_axis,
                                     y_axis,
                                     global_across_channels=TRUE,
                                     palette=c("blue","green","red"),
                                     resolution=72,
                                     data_transformation_reverse=identity,
                                     file_name="biplot.pdf",
                                     raster.height=480,
                                     raster.width=480,
                                     ... #pass to plot for e.g. tSNE biplot
                                     )
{
    regular_channels <- setdiff(colnames(matrix),c(x_axis,y_axis))

    if(global_across_channels){
        data_range <- range(matrix[,regular_channels],na.rm=TRUE)
        data_range <- matrix(rep(data_range,length(regular_channels)),ncol=length(regular_channels),nrow=2,byrow=FALSE,dimnames=list(c("min","max"),regular_channels))
    } else {
        data_range <- apply(matrix[,regular_channels],na.rm=TRUE,2,range,na.rm=TRUE)
        rownames(data_range) <- c("min","max")
    }

    x <- matrix[,x_axis]
    y <- matrix[,y_axis]
    xp <- matrix[,regular_channels,drop=FALSE]

    if(any(!(is.na(x)&is.na(y)))){
        rasters <- lapply(
            regular_channels,
            function(pname,xp,data.range,x,y)
            {
                color.scale <- unique(colorRampPalette(palette)(1000))
                n <- length(color.scale)

                breaks <- unique(seq(data.range["min",pname],data.range["max",pname],length.out=n+1))
                if(length(unique(breaks))>1){
                    points.colors <- as.character(cut(xp[,pname],breaks=breaks,labels=color.scale))
                } else {
                    points.colors <- rep("lightsteelblue",length(xp[,pname]))
                }
                mainplot <- paste(tmpDir(),"/mainplot_",pname,".png",sep="")
                png(mainplot,res=resolution,height=raster.height*resolution/72,width=raster.width*resolution/72)
                par("bty"="l")
                plot(
                    x,
                    y,
                    col=points.colors,
                    xlab=x_axis,
                    ylab=y_axis,
                    main=pname,
                    ...
                )
                dev.off()

                colorscale <- paste(tmpDir(),"/colorscale_",pname,".png",sep="")
                png(colorscale,res=resolution,height=raster.height/2*resolution/72,width=raster.width*resolution/72)
                plot.new()
                par("mar"=c(2,1,2,1))


                xlims <- par("usr")[c(1, 2)]
                x_coords <- seq(xlims[1],xlims[2],length.out=n+1)
                ylims <- par("usr")[c(3, 4)]

                labels <- signif(data_transformation_reverse(breaks),2)
                labels <- labels[round(seq(1,length(labels),length.out=5))]
                labels.x_coords <- seq(x_coords[1],x_coords[length(x_coords)],length.out=5)

                rect(border=NA,ybottom=ylims[1],ytop=ylims[2],xleft=x_coords[-length(x_coords)],xright=x_coords[-1],col=color.scale)
                text(xpd=TRUE,y=ylims[1],pos=1,labels=labels,x=labels.x_coords)
                text(xpd=TRUE,y=ylims[2],pos=3,labels=paste(pname,"intensity"),x=mean(xlims))
                dev.off()

                return(list(main.file=mainplot,scale.file=colorscale))
            },
            xp=xp,
            data.range=data_range,
            x=x,
            y=y
        )
        
        file <- file_name

        pdf(file)
        if(global_across_channels){
            plot.new()
            grid.raster(readPNG(rasters[[1]]$scale.file,native=TRUE))
        }
        lapply(
            rasters,
            function(x){
                par("mar"=c(0,0,0,0))
                grid.newpage()
                label <- sub(".png","",sub("mainplot_","",tail(strsplit(x$main.file,"/")[[1]],1),fixed=TRUE),fixed=TRUE)
                
                if(!global_across_channels){
                    grid.raster(readPNG(x$main.file,native=TRUE),y=0.6,height=0.8)
                    grid.raster(readPNG(x$scale.file,native=TRUE),y=0.1,height=0.2)
                }
                if(global_across_channels){
                    grid.raster(readPNG(x$main.file,native=TRUE))
                }
                grid.text(x=unit(1,"npc"),y=unit(1,"npc"),label=label,just=c(1,1),gp=gpar(col="white",cex=0.1))
                return(NULL)
            }
        )
        dev.off()
    }
    NULL
}

utils::globalVariables(c("predictions_leave_one_out_backbone", "predictions_variable_targets", "events.code", "predictors")) ## These are passed from the predictions object and assigned from its named list return
