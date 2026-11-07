close all

load("multi_antenna_base_results.mat");

params_latex = ["\hat{\alpha}","\hat{t}_0"];
params_text = ["alpha","t0"];

RGB = orderedcolors("gem12");
H = compose("#%02X%02X%02X",round(RGB*255));
colors = cellstr(H);
colorMap = containers.Map(cellstr([unique(iter_arr)]), colors(1:numel([unique(iter_arr)])));

C = orderedcolors("gem12");

gem12 = containers.Map( ...
    {'blue', 'orange', 'yellow', 'purple', 'green', ...
    'cyan', 'magenta', 'darkBlue', 'darkOrange', ...
    'darkYellow', 'darkPurple', 'darkGreen'}, ...
    num2cell(C, 2));

exclude = all_strats(contains(all_strats, "hetero.") & contains(all_strats, "prop. scaling")); % ["None"];

% Conventional baselines in this results file (set by conventional_methods in multi_antenna.m)
conv_names = string(fieldnames(avg_dep_mse_conv)).';
conv_style = containers.Map({'median', 'trimmed_mean', 'skipped_mean', 'skipped_wmean'}, ...
    {{gem12("orange"), ':x', 'median'}, {gem12("darkGreen"), ':^', 'trimmed mean'}, ...
     {gem12("cyan"), ':s', 'skipped mean'}, {gem12("darkPurple"), ':d', 'weighted skipped mean'}});

for agent_db_idx = 1:length(agent_db_values)
    for atk_idx = 1:numel(attacks)
        atk_type = attacks(atk_idx);
        selected_atk_idxs = 1:numel(attacks);
        selected_atks = attacks;
        for param_idx = 1:2
            % normalize performance metrics
            if param_idx == 1
                norm_coeff = alpha_true^2;
            else
                norm_coeff = t0_true^2;
            end
    
            avg_dep_crlb(1,agent_db_idx,:,param_idx,:,atk_idx) = avg_dep_crlb(1,agent_db_idx,:,param_idx,:,atk_idx) ./ norm_coeff;
            avg_dep_var(1,agent_db_idx,:,param_idx,:,atk_idx) = avg_dep_var(1,agent_db_idx,:,param_idx,:,atk_idx) ./ norm_coeff;
            avg_dep_mse(1,agent_db_idx,:,param_idx,:,atk_idx) = avg_dep_mse(1,agent_db_idx,:,param_idx,:,atk_idx) ./ norm_coeff;
            if attacker_enabled
                avg_dep_mse_undefended(1,agent_db_idx,:,param_idx,:,atk_idx) = avg_dep_mse_undefended(1,agent_db_idx,:,param_idx,:,atk_idx) ./ norm_coeff;
                avg_dep_mse_oracle(1,agent_db_idx,:,param_idx,:,atk_idx) = avg_dep_mse_oracle(1,agent_db_idx,:,param_idx,:,atk_idx) ./ norm_coeff;
                avg_dep_mse_full_honest(1,agent_db_idx,:,param_idx,:,atk_idx) = avg_dep_mse_full_honest(1,agent_db_idx,:,param_idx,:,atk_idx) ./ norm_coeff;
                for cn = conv_names
                    avg_dep_mse_conv.(char(cn))(1,agent_db_idx,:,param_idx,:,atk_idx) = avg_dep_mse_conv.(char(cn))(1,agent_db_idx,:,param_idx,:,atk_idx) ./ norm_coeff;
                end
                % avg_dep_mse_basis_supp(1,agent_db_idx,:,param_idx,:,atk_idx) = avg_dep_mse_basis_supp(1,agent_db_idx,:,param_idx,:,atk_idx) ./ norm_coeff;
                avg_dep_crlb_oracle(1,agent_db_idx,:,param_idx,:,atk_idx) = avg_dep_crlb_oracle(1,agent_db_idx,:,param_idx,:,atk_idx) ./ norm_coeff;
            end
    
            % Collect all plot values for y-axis scaling.
            all_vals = [avg_dep_mse(1,agent_db_idx,:,param_idx,:,atk_idx);
                        % avg_dep_crlb(selected_strat_idxs,agent_db_idx,:,param_idx,:,attack_idx);
                        % avg_dep_var(selected_strat_idxs,agent_db_idx,:,param_idx,:,attack_idx)
                        ];
            if attacker_enabled
                all_vals = [all_vals; avg_dep_mse_undefended(1,agent_db_idx,:,param_idx,:,atk_idx);
                            avg_dep_mse_oracle(1,agent_db_idx,:,param_idx,:,atk_idx);
                            avg_dep_mse_full_honest(1,agent_db_idx,:,param_idx,:,atk_idx)];
                for cn = conv_names
                    all_vals = [all_vals; avg_dep_mse_conv.(char(cn))(1,agent_db_idx,:,param_idx,:,atk_idx)];
                end
                            % avg_dep_mse_basis_supp(1,agent_db_idx,:,param_idx,:,atk_idx)];
                            % avg_dep_crlb_oracle(selected_strat_idxs,agent_db_idx,:,param_idx,:,attack_idx)];
            end
    
            for scheme_idx = 1:length(selected_schemes)
                scheme = selected_schemes(scheme_idx);
    
                count_vals = sensor_vals;
    
                for constraint_idx = 1:length(constraints)
                    constraint = constraints(constraint_idx);
                    data_idx_seq = selected_atk_idxs; % selected_atk_idxs(find(contains(selected_atks, constraint)).');

                    fig1 = figure("Name", atk_type);
                    plot_w = 450;
                    set(fig1,'Position',[100,100,plot_w,plot_w])
                    ax = gca;
                    hold(ax, 'on');

                    x_axis_series = channel_db_values(scheme_idx,:,agent_db_idx);
                    plot_line_width = 1.5;
    
                    plot(ax, x_axis_series,(squeeze(avg_dep_mse(1,agent_db_idx,:,param_idx,scheme_idx,atk_idx))),'-x','color',gem12("blue"),'LineWidth', plot_line_width, 'DisplayName', 'proposed')
                    
                    % plot(ax, x_axis_series,(squeeze(avg_dep_var(strat_idx,agent_db_idx,:,param_idx,scheme_idx,atk_idx,1))),'-^','color',colorMap(iter_key),'LineWidth', plot_line_width)
                    % plot(ax, x_axis_series,(squeeze(avg_dep_crlb(strat_idx,agent_db_idx,:,param_idx,scheme_idx,atk_idx,1))),'--o','color',colorMap(iter_key),'LineWidth', plot_line_width)
                    if attacker_enabled
                        plot(ax, x_axis_series,(squeeze(avg_dep_mse_undefended(1,agent_db_idx,:,param_idx,scheme_idx,atk_idx))),':x','color',gem12("purple"),'LineWidth', plot_line_width, 'DisplayName', 'undefended')
                        plot(ax, x_axis_series,(squeeze(avg_dep_mse_full_honest(1,agent_db_idx,:,param_idx,scheme_idx,atk_idx))),':x','color',gem12("magenta"),'LineWidth', plot_line_width, 'DisplayName', 'attack-free')
                        for cn = conv_names
                            st = conv_style(char(cn));   % {colour, line style, legend name}
                            plot(ax, x_axis_series,(squeeze(avg_dep_mse_conv.(char(cn))(1,agent_db_idx,:,param_idx,scheme_idx,atk_idx))), ...
                                st{2},'color',st{1},'LineWidth', plot_line_width, 'DisplayName', st{3})
                        end
                        plot(ax, x_axis_series,(squeeze(avg_dep_mse_oracle(1,agent_db_idx,:,param_idx,scheme_idx,atk_idx))),':x','color',gem12("green"),'LineWidth', plot_line_width, 'DisplayName', 'oracle')
                        
                        % plot(ax, x_axis_series,(squeeze(avg_dep_mse_basis_supp(1,agent_db_idx,:,param_idx,scheme_idx,atk_idx,1))),':x','color',gem12("green"),'LineWidth', plot_line_width)
                        % plot(ax, x_axis_series,(squeeze(avg_dep_crlb_oracle(1,agent_db_idx,:,param_idx,scheme_idx,atk_idx,1))),':o','color','magenta','LineWidth', plot_line_width)
                    end
    
                    xlabel('Channel SNR (dB)')
                    ylabel(" ")
    
                    title('$$'+params_latex(param_idx)+'$$','Interpreter','latex')
    
                    ax.FontSize = 15;          % font size
                    ax.YScale = 'log';          % log scale if needed
                    grid(ax, 'on');             % turn on grid
                    box(ax, 'on');              % border on all four sides

                    y_lower = 10^(-0.1)*min(all_vals,[],"all");
                    y_upper = 1.05*max(all_vals,[],"all");
                    ylim([y_lower,y_upper]);
    
                    ax.XTick = channel_snr;

                    if param_idx == 1
                        legend("Location","best")
                    end

                    home_dir = char(java.lang.System.getProperty('user.home'));
                    folder_path = home_dir + "\Desktop\Byzantine AirComp Results\";

                    save_folder = folder_path + experiment;
                    if ~exist(save_folder, 'dir')
                        mkdir(save_folder);
                    end
                    exportgraphics(fig1, fullfile(save_folder, sprintf(experiment + "_" + atk_type + "_" + string(agent_db_values(agent_db_idx)) + "_db_Tp=" + string(Tp) + "_" + params_text(param_idx) + "_" + selected_schemes(scheme_idx) + ".png")), 'Resolution', 300);
                end
            end
        end
    end
end