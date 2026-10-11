% close all

load("multi_antenna_disjoint_results.mat");

params_latex = ["\hat{\alpha}","\hat{t}_0"];
params_text = ["alpha","t0"];

RGB = orderedcolors("gem12");
H = compose("#%02X%02X%02X",round(RGB*255));
colors = cellstr(H);
colorMap = containers.Map(cellstr([unique(iter_arr)]), colors(1:numel([unique(iter_arr)])));

exclude = all_strats(contains(all_strats, "hetero.") & contains(all_strats, "prop. scaling")); % ["None"];

for atk_idx = 1:numel(attacks)
    atk_type = attacks(atk_idx);
    selected_atk_idxs = 1;
    selected_atks = attacks(selected_atk_idxs);
    for agent_db_idx = 1:length(agent_db_values)
        for param_idx = 1:2
            % normalize performance metrics
            if param_idx == 1
                norm_coeff = alpha_true^2;
            else
                norm_coeff = t0_true^2;
            end
    
            avg_dep_crlb(selected_atk_idxs,agent_db_idx,:,param_idx,:,atk_idx) = avg_dep_crlb(selected_atk_idxs,agent_db_idx,:,param_idx,:,atk_idx) ./ norm_coeff;
            avg_dep_var(selected_atk_idxs,agent_db_idx,:,param_idx,:,atk_idx) = avg_dep_var(selected_atk_idxs,agent_db_idx,:,param_idx,:,atk_idx) ./ norm_coeff;
            avg_dep_mse(selected_atk_idxs,agent_db_idx,:,param_idx,:,atk_idx) = avg_dep_mse(selected_atk_idxs,agent_db_idx,:,param_idx,:,atk_idx) ./ norm_coeff;
            if attacker_enabled
                avg_dep_mse_undefended(selected_atk_idxs,agent_db_idx,:,param_idx,:,atk_idx) = avg_dep_mse_undefended(selected_atk_idxs,agent_db_idx,:,param_idx,:,atk_idx) ./ norm_coeff;
                avg_dep_mse_oracle(selected_atk_idxs,agent_db_idx,:,param_idx,:,atk_idx) = avg_dep_mse_oracle(selected_atk_idxs,agent_db_idx,:,param_idx,:,atk_idx) ./ norm_coeff;
                avg_dep_mse_full_honest(selected_atk_idxs,agent_db_idx,:,param_idx,:,atk_idx) = avg_dep_mse_full_honest(selected_atk_idxs,agent_db_idx,:,param_idx,:,atk_idx) ./ norm_coeff;
                avg_dep_mse_trimmed(selected_atk_idxs,agent_db_idx,:,param_idx,:,atk_idx) = avg_dep_mse_trimmed(selected_atk_idxs,agent_db_idx,:,param_idx,:,atk_idx) ./ norm_coeff;
                % avg_dep_mse_basis_supp(selected_atk_idxs,agent_db_idx,:,param_idx,:,atk_idx) = avg_dep_mse_basis_supp(selected_atk_idxs,agent_db_idx,:,param_idx,:,atk_idx) ./ norm_coeff;
                avg_dep_crlb_oracle(selected_atk_idxs,agent_db_idx,:,param_idx,:,atk_idx) = avg_dep_crlb_oracle(selected_atk_idxs,agent_db_idx,:,param_idx,:,atk_idx) ./ norm_coeff;
            end
    
            % Collect all plot values for y-axis scaling.
            all_vals = [avg_dep_mse(selected_atk_idxs,agent_db_idx,:,param_idx,:,atk_idx);
                        % avg_dep_crlb(selected_strat_idxs,agent_db_idx,:,param_idx,:,attack_idx);
                        % avg_dep_var(selected_strat_idxs,agent_db_idx,:,param_idx,:,attack_idx)
                        ];
            if attacker_enabled
                all_vals = [all_vals; avg_dep_mse_undefended(selected_atk_idxs,agent_db_idx,:,param_idx,:,atk_idx);
                            avg_dep_mse_oracle(selected_atk_idxs,agent_db_idx,:,param_idx,:,atk_idx);
                            avg_dep_mse_full_honest(selected_atk_idxs,agent_db_idx,:,param_idx,:,atk_idx);
                            avg_dep_mse_trimmed(selected_atk_idxs,agent_db_idx,:,param_idx,:,atk_idx)];
                            % avg_dep_mse_basis_supp(selected_atk_idxs,agent_db_idx,:,param_idx,:,atk_idx)];
                            % avg_dep_crlb_oracle(selected_strat_idxs,agent_db_idx,:,param_idx,:,attack_idx)];
            end
    
            for scheme_idx = 1:length(selected_schemes)
                scheme = selected_schemes(scheme_idx);
    
                count_vals = sensor_vals;
    
                for constraint_idx = 1:length(constraints)
                    constraint = constraints(constraint_idx);
                    data_idx_seq = selected_atk_idxs; % selected_atk_idxs(find(contains(selected_atks, constraint)).');
    
                    % if attack_idx == 1
                    %     figLgd = figure;
                    %     axLgd = gca;
                    %     plot_w = 450;
                    %     set(figLgd,'Position',[100,100,plot_w,plot_w])
                    %     hold(axLgd, 'on');
                    % 
                    %     lgd = [];
                    % 
                    %     lgd = ["proposed", "oracle", "undefended", "full honest", "trimmed"];
                    % 
                    %     plot(axLgd,nan,nan,'color',gem12("blue"),'LineWidth',2);
                    %     plot(axLgd,nan,nan,'color','green','LineWidth',2);
                    %     plot(axLgd,nan,nan,'color',gem12("purple"),'LineWidth',2);
                    %     plot(axLgd,nan,nan,'color',gem12("magenta"),'LineWidth',2);
                    %     plot(axLgd,nan,nan,'color',gem12("orange"),'LineWidth',2);
                    %     % plot(axLgd,nan,nan,'color',gem12("green"),'LineWidth',2);
                    % 
                    %     plot(axLgd, nan,nan,'x','color','black');
                    %     % plot(axLgd, nan,nan,'^','color','black');
                    %     % plot(axLgd, nan,nan,'o','color','black');
                    % 
                    %     lgdObj = legend(axLgd, [lgd, "MSE"]);
                    %     lgdObj.Location = 'none';
                    %     lgdObj.Units = 'normalized';
                    %     lgdObj.Position(1) = (1 - lgdObj.Position(3)) / 2;  % center horizontally
                    %     lgdObj.Position(2) = (1 - lgdObj.Position(4)) / 2;  % center vertically
                    % 
                    %     lgdObj.FontSize = 15;        % increase text size
                    %     lgdObj.ItemTokenSize = [30 18];  % increase marker/line size [width height]
                    %     axLgd.Color = 'none';
                    %     axLgd.XColor = 'none';
                    %     axLgd.YColor = 'none';
                    %     figLgd.Color = 'white';
                    %     drawnow;
                    % end

                    fig1 = figure("Name", atk_type);
                    plot_w = 450;
                    set(fig1,'Position',[100,100,plot_w,plot_w])
                    ax = gca;
                    hold(ax, 'on');

                    x_axis_series = channel_db_values(scheme_idx,:,agent_db_idx);
                    plot_line_width = 1.5;
    
                    for count_idx = 1:length(count_vals)
                        for strat_idx = data_idx_seq
                            strat = all_strats(strat_idx);
                            split_str = split(strat);
                            iter_key = join(split_str(2:end), " ");
    
                            plot(ax, x_axis_series,(squeeze(avg_dep_mse(strat_idx,agent_db_idx,:,param_idx,scheme_idx,count_idx,1))),'-x','color',gem12("blue"),'LineWidth', plot_line_width)
                            % plot(ax, x_axis_series,(squeeze(avg_dep_var(strat_idx,agent_db_idx,:,param_idx,scheme_idx,count_idx,1))),'-^','color',colorMap(iter_key),'LineWidth', plot_line_width)
                            % plot(ax, x_axis_series,(squeeze(avg_dep_crlb(strat_idx,agent_db_idx,:,param_idx,scheme_idx,count_idx,1))),'--o','color',colorMap(iter_key),'LineWidth', plot_line_width)
                            if attacker_enabled
                                plot(ax, x_axis_series,(squeeze(avg_dep_mse_undefended(strat_idx,agent_db_idx,:,param_idx,scheme_idx,count_idx,1))),':x','color',gem12("purple"),'LineWidth', plot_line_width)
                                plot(ax, x_axis_series,(squeeze(avg_dep_mse_full_honest(strat_idx,agent_db_idx,:,param_idx,scheme_idx,count_idx,1))),':x','color',gem12("magenta"),'LineWidth', plot_line_width)
                                plot(ax, x_axis_series,(squeeze(avg_dep_mse_trimmed(strat_idx,agent_db_idx,:,param_idx,scheme_idx,count_idx,1))),':x','color',gem12("orange"),'LineWidth', plot_line_width)
                                % plot(ax, x_axis_series,(squeeze(avg_dep_mse_basis_supp(strat_idx,agent_db_idx,:,param_idx,scheme_idx,count_idx,1))),':x','color',gem12("green"),'LineWidth', plot_line_width)
                                plot(ax, x_axis_series,(squeeze(avg_dep_mse_oracle(strat_idx,agent_db_idx,:,param_idx,scheme_idx,count_idx,1))),':x','color','green','LineWidth', plot_line_width)
                                % plot(ax, x_axis_series,(squeeze(avg_dep_crlb_oracle(strat_idx,agent_db_idx,:,param_idx,scheme_idx,count_idx,1))),':o','color','magenta','LineWidth', plot_line_width)
                            end
                        end
                    end
    
                    xlabel('Channel SNR (dB)')
                    ylabel(" ")
    
                    title('$$'+params_latex(param_idx)+'$$','Interpreter','latex')
    
                    ax.FontSize = 15;          % font size
                    ax.YScale = 'log';          % log scale if needed
                    grid(ax, 'on');             % turn on grid
    
                    y_lower = 10^(-0.1)*min(all_vals,[],"all");
                    y_upper = 1.05*max(all_vals,[],"all");
                    ylim([y_lower,y_upper]);
    
                    ax.XTick = channel_snr;
                end
            end
        end
    end
end