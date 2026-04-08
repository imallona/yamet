"""
Coverage and variability simulations.

Source code (Rmd scripts, R helpers) lives in workflow/rules/src/ with a
simulations_ prefix. Outputs go under workflow/simulations/ (output/ for yamet
data, results/ for HTML reports, schemes/ for figures, rds/ for intermediates).

The gastrulation data needed by the coverage simulations is downloaded once by
the download_argelaguet rule in argelaguet.smk. Simulation rules depend on
argelaguet/downloaded.flag directly, with no symlinks.
"""

import csv

SIM_BASE = "simulations"
SIM_INPUT = op.join(SIM_BASE, "input")
SIM_OUTPUT = op.join(SIM_BASE, "output")
SIM_RESULTS = op.join(SIM_BASE, "results")
SIM_SCHEMES = op.join(SIM_BASE, "schemes")
SIM_RDS = op.join(SIM_BASE, "rds")

_SIM_YAMET = op.join(op.expanduser("~"), ".local", "yamet", "bin", "yamet")

WC_S_GRID = {"S": [20, 40, 60, 80, 100], "N": 300, "f": 641}
WC_N_GRID = {"S": 30, "N": [100, 200, 300, 400, 500], "f": 641}
WC_DISP = {"S": 20, "N": 500, "f": 321}

BC_S_GRID = {"S": [20, 40, 60, 80, 100], "N": 200, "f": 100}
BC_N_GRID = {"S": 30, "N": [100, 200, 300, 400, 500], "f": 100}
BC_DISP = {"S": 20, "N": 500, "f": 200}


def _write_sim_interval_bed(parameters_path, file_path, chrom="chrSim", start=1):
    df = pd.read_table(parameters_path, delimiter="\t", header=0)
    end = df["n_cpgs"].max()
    with open(file_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([chrom, str(start), str(end)])


## -- emanuel's coverage simulations ------------------------------------------

rule sim_simulate_data:
    input:
        parameters_path=op.join(SIM_INPUT, "parameters.tsv"),
        downloaded_data=op.join(ARGELAGUET_BASE, "downloaded.flag"),
    output:
        directory(op.join(SIM_OUTPUT, "sim_data")),
        op.join(SIM_OUTPUT, "sim_data",
                "sim_cell_9_50_50_rand_medium_lmr.tsv"),
        html=op.join(SIM_RESULTS, "simulation_01_sim_data.html"),
    log:
        op.join("logs", "sim_simulate_data.log"),
    conda:
        op.join("..", "envs", "sim_r_emanuel.yaml")
    params:
        rmd=op.abspath(op.join("rules", "src", "simulations_01_sim_data.Rmd")),
        parameters_abs=lambda wc: op.abspath(op.join(SIM_INPUT, "parameters.tsv")),
        low_real_dir=lambda wc: op.abspath(
            op.join(ARGELAGUET_BASE, "met", "cpg_level")),
        out_dir=lambda wc: op.abspath(op.join(SIM_OUTPUT, "sim_data")),
        results_abs=lambda wc: op.abspath(SIM_RESULTS),
    shell:
        """
        mkdir -p {SIM_OUTPUT} {SIM_RESULTS}
        Rscript -e 'rmarkdown::render("{params.rmd}",
                        "html_document",
                        output_file=file.path("{params.results_abs}", "simulation_01_sim_data.html"),
                        params=list(
                            parameters_path="{params.parameters_abs}",
                            low_real_dir="{params.low_real_dir}",
                            out_dir="{params.out_dir}"))' \
            &> {log}
        """


rule sim_visualize_coverage:
    input:
        sim_data_file=op.join(SIM_RESULTS, "simulation_01_sim_data.html"),
        sim_data_dir=op.join(SIM_OUTPUT, "sim_data"),
        flag=op.join(SIM_OUTPUT, "sim_data",
                     "sim_cell_9_50_50_rand_medium_lmr.tsv"),
    output:
        html=op.join(SIM_RESULTS, "simulation_02_coverage_plots.html"),
        coverage_data=op.join(SIM_OUTPUT, "sim_data_coverage",
                              "overall_coverage.rds"),
        stretches_length=op.join(SIM_OUTPUT, "sim_data_coverage",
                                 "covered_stretches_length.rds"),
    log:
        op.join("logs", "sim_visualize_coverage.log"),
    conda:
        op.join("..", "envs", "sim_r_plot.yaml")
    params:
        rmd=op.abspath(op.join("rules", "src", "simulations_02_coverage_plots.Rmd")),
        sim_data_dir=lambda wc: op.abspath(op.join(SIM_OUTPUT, "sim_data")),
        coverage_out_dir=lambda wc: op.abspath(op.join(SIM_OUTPUT, "sim_data_coverage")),
        fig_path=lambda wc: op.abspath(op.join(SIM_OUTPUT, "02_plots")) + "/",
        results_abs=lambda wc: op.abspath(SIM_RESULTS),
    shell:
        """
        Rscript -e 'rmarkdown::render("{params.rmd}",
                        "html_document",
                        output_file=file.path("{params.results_abs}", "simulation_02_coverage_plots.html"),
                        params=list(
                            sim_data_dir="{params.sim_data_dir}",
                            coverage_out_dir="{params.coverage_out_dir}",
                            fig_path="{params.fig_path}"))' \
            &> {log}
        """


rule sim_write_reference_file:
    input:
        parameters_file=op.join(SIM_INPUT, "parameters.tsv"),
        sim_data=op.join(SIM_OUTPUT, "sim_data"),
        sim_data_file=op.join(SIM_RESULTS, "simulation_01_sim_data.html"),
    output:
        bed=op.join(SIM_BASE, "intervals.bed"),
    log:
        op.join("logs", "sim_write_reference_file.log"),
    run:
        _write_sim_interval_bed(
            parameters_path=input.parameters_file,
            file_path=output.bed,
        )


rule sim_run_yamet:
    conda:
        op.join("..", "envs", "sim_yamet.yaml")
    input:
        yamet=_SIM_YAMET,
        sim_data=op.join(SIM_OUTPUT, "sim_data"),
        flag=op.join(SIM_OUTPUT, "sim_data",
                     "sim_cell_9_50_50_rand_medium_lmr.tsv"),
        sim_data_file=op.join(SIM_RESULTS, "simulation_01_sim_data.html"),
        bed=op.join(SIM_BASE, "intervals.bed"),
    output:
        out=directory(op.join(SIM_OUTPUT, "yamet", "out")),
        flag=op.join(SIM_OUTPUT, "yamet", "out", "10000.tsv"),
    log:
        op.join("logs", "sim_run_yamet.log"),
    threads: 6
    params:
        ncpgs="50 100 1000 10000",
    shell:
        """
        mkdir -p {output.out}
        for ncpg in {params.ncpgs};
        do
           {input.yamet} \
               --cell {input.sim_data}/sim_cell_*_*_${{ncpg}}_*.tsv \
               --reference {input.sim_data}/cpgPositions_${{ncpg}}.tsv \
               --intervals {input.bed} \
               --out {output.out}/${{ncpg}}.tsv \
               --norm-det-out {output.out}/${{ncpg}}.det_norm \
               --det-out {output.out}/${{ncpg}}.det \
               --cores {threads} \
               --no-print-sampens
        done &> {log}
        """


rule sim_visualize_yamet_results:
    conda:
        op.join("..", "envs", "sim_r_plot.yaml")
    input:
        covplots=op.join(SIM_RESULTS, "simulation_02_coverage_plots.html"),
        yamet_dir=op.join(SIM_OUTPUT, "yamet", "out"),
        yamet=_SIM_YAMET,
        flag=op.join(SIM_OUTPUT, "sim_data",
                     "sim_cell_9_50_50_rand_medium_lmr.tsv"),
        another_flag=op.join(SIM_OUTPUT, "yamet", "out", "10000.tsv"),
    output:
        html=op.join(SIM_RESULTS, "simulation_03_yamet_plots.html"),
    log:
        op.join("logs", "sim_visualize_yamet_results.log"),
    params:
        rmd=op.abspath(op.join("rules", "src", "simulations_03_yamet_plots.Rmd")),
        yamet_dir=lambda wc: op.abspath(op.join(SIM_OUTPUT, "yamet", "out")),
        fig_path=lambda wc: op.abspath(op.join(SIM_OUTPUT, "03_plots")) + "/",
        results_abs=lambda wc: op.abspath(SIM_RESULTS),
    shell:
        """
        Rscript -e 'rmarkdown::render("{params.rmd}",
                        "html_document",
                        output_file=file.path("{params.results_abs}", "simulation_03_yamet_plots.html"),
                        params=list(
                            yamet_dir="{params.yamet_dir}",
                            fig_path="{params.fig_path}"))' \
            &> {log}
        """


rule sim_compile_figure_2:
    conda:
        op.join("..", "envs", "sim_r_plot.yaml")
    input:
        yamet=_SIM_YAMET,
        yamet_dir=op.join(SIM_OUTPUT, "yamet", "out"),
        yamet_plots=op.join(SIM_RESULTS, "simulation_03_yamet_plots.html"),
        coverage_data=op.join(SIM_OUTPUT, "sim_data_coverage",
                              "overall_coverage.rds"),
        stretches_length=op.join(SIM_OUTPUT, "sim_data_coverage",
                                 "covered_stretches_length.rds"),
    output:
        fig_pdf=op.join(SIM_SCHEMES, "Figure2.pdf"),
        fig_svg=op.join(SIM_SCHEMES, "Figure2.svg"),
        html=op.join(SIM_RESULTS, "simulation_figure2.html"),
    log:
        op.join("logs", "sim_compile_figure_2.log"),
    params:
        rmd=op.abspath(op.join("rules", "src", "simulations_figure2.Rmd")),
        yamet_dir=lambda wc: op.abspath(op.join(SIM_OUTPUT, "yamet", "out")),
        schemes_dir=lambda wc: op.abspath(SIM_SCHEMES),
        coverage_data_dir=lambda wc: op.abspath(op.join(SIM_OUTPUT, "sim_data_coverage")),
        results_abs=lambda wc: op.abspath(SIM_RESULTS),
    shell:
        """
        mkdir -p {SIM_SCHEMES}
        Rscript -e 'rmarkdown::render("{params.rmd}",
                        "html_document",
                        output_file=file.path("{params.results_abs}", "simulation_figure2.html"),
                        params=list(
                            yamet_dir="{params.yamet_dir}",
                            schemes_dir="{params.schemes_dir}",
                            coverage_data_dir="{params.coverage_data_dir}"))' \
            &> {log}
        """


## -- feature-level differential simulations -----------------------------------

rule sim_feat_diff_template:
    conda:
        op.join("..", "envs", "sim_r_sim.yaml")
    output:
        op.join(SIM_OUTPUT, "feat_diff", "intervals.{N}.{f}.tsv"),
    log:
        op.join("logs", "sim_feat_diff_template_{N}_{f}.log"),
    script:
        op.join("src", "simulations_feat_diff_createTemplate.R")


rule sim_feat_diff_samples:
    conda:
        op.join("..", "envs", "sim_r_sim.yaml")
    input:
        op.join(SIM_OUTPUT, "feat_diff", "intervals.{N}.{f}.tsv"),
    output:
        touch(op.join(SIM_OUTPUT, "feat_diff", "{N}.{f}.flag")),
    log:
        op.join("logs", "sim_feat_diff_samples_{N}_{f}.log"),
    params:
        samples=50,
        data=op.join(SIM_OUTPUT, "feat_diff"),
    script:
        op.join("src", "simulations_feat_diff_getSamples.R")


rule sim_feat_diff_yamet:
    conda:
        op.join("..", "envs", "sim_yamet.yaml")
    input:
        op.join(SIM_OUTPUT, "feat_diff", "{N}.{f}.flag"),
        yamet=_SIM_YAMET,
    output:
        det_out=op.join(SIM_OUTPUT, "feat_diff",
                        "yamet.{N}.{f}.det.out"),
        norm_det_out=op.join(SIM_OUTPUT, "feat_diff",
                             "yamet.{N}.{f}.norm.det.out"),
    log:
        op.join("logs", "sim_feat_diff_yamet_{N}_{f}.log"),
    params:
        data=op.join(SIM_OUTPUT, "feat_diff"),
    threads: 3
    shell:
        """
        {input.yamet} \
            --cell {params.data}/sim.*.tsv \
            --reference {params.data}/ref.{wildcards.N}.{wildcards.f}.tsv \
            --intervals {params.data}/intervals.{wildcards.N}.{wildcards.f}.tsv \
            --det-out {output.det_out} \
            --norm-det-out {output.norm_det_out} \
            --meth-out {params.data}/yamet.{wildcards.N}.{wildcards.f}.meth.out \
            --cores {threads} \
            --no-print-sampens \
            &> {log}
        """


rule sim_feat_diff_plots:
    conda:
        op.join("..", "envs", "sim_r_plot.yaml")
    input:
        det_out=op.join(SIM_OUTPUT, "feat_diff",
                        "yamet.500.200.det.out"),
        norm_det_out=op.join(SIM_OUTPUT, "feat_diff",
                             "yamet.500.200.norm.det.out"),
    params:
        N=lambda wildcards, input: op.basename(input.det_out).split(".")[1],
        f=lambda wildcards, input: op.basename(input.det_out).split(".")[2],
    output:
        html=op.join(SIM_RESULTS, "simulation_04_feat_diff.html"),
    log:
        op.join("logs", "sim_feat_diff_plots.log"),
    script:
        op.join("src", "simulations_04_feat_diff.Rmd")


## -- within-cell variability simulations --------------------------------------

rule sim_within_cell_template:
    conda:
        op.join("..", "envs", "sim_r_sim.yaml")
    output:
        op.join(SIM_OUTPUT, "within_cell",
                "intervals.{S}.{N}.{f}.tsv"),
    log:
        op.join("logs", "sim_within_cell_template_{S}_{N}_{f}.log"),
    script:
        op.join("src", "simulations_within_cell_createTemplate.R")


rule sim_within_cell_samples:
    conda:
        op.join("..", "envs", "sim_r_sim.yaml")
    input:
        op.join(SIM_OUTPUT, "within_cell",
                "intervals.{S}.{N}.{f}.tsv"),
    output:
        op.join(SIM_OUTPUT, "within_cell",
                "sim.{s}.{S}.{N}.{f}.tsv"),
    log:
        op.join("logs", "sim_within_cell_samples_{s}_{S}_{N}_{f}.log"),
    params:
        data=op.join(SIM_OUTPUT, "within_cell"),
    script:
        op.join("src", "simulations_within_cell_getSamples.R")


rule sim_within_cell_yamet:
    conda:
        op.join("..", "envs", "sim_yamet.yaml")
    input:
        cells=lambda wildcards: expand(
            op.join(SIM_OUTPUT, "within_cell",
                    "sim.{s}.{{S}}.{{N}}.{{f}}.tsv"),
            s=range(1, int(wildcards.S) + 1),
        ),
        yamet=_SIM_YAMET,
    output:
        det_out=op.join(SIM_OUTPUT, "within_cell",
                        "yamet.{S}.{N}.{f}.det.out"),
        norm_det_out=op.join(SIM_OUTPUT, "within_cell",
                             "yamet.{S}.{N}.{f}.norm.det.out"),
    log:
        op.join("logs", "sim_within_cell_yamet_{S}_{N}_{f}.log"),
    params:
        data=op.join(SIM_OUTPUT, "within_cell"),
    threads: 4
    benchmark:
        repeat(op.join(SIM_OUTPUT, "within_cell",
                       "yamet.{S}.{N}.{f}.benchmark"), 5)
    shell:
        """
        {input.yamet} \
            --cell {input.cells} \
            --reference {params.data}/ref.{wildcards.N}.{wildcards.f}.tsv \
            --intervals {params.data}/intervals.{wildcards.S}.{wildcards.N}.{wildcards.f}.tsv \
            --skip-header-intervals \
            --det-out {output.det_out} \
            --norm-det-out {output.norm_det_out} \
            --meth-out {params.data}/yamet.{wildcards.S}.{wildcards.N}.{wildcards.f}.meth.out \
            --cores {threads} \
            --no-print-sampens \
            &> {log}
        """


rule sim_within_cell_scmet:
    conda:
        op.join("..", "envs", "sim_r_plot.yaml")
    input:
        op.join(SIM_OUTPUT, "within_cell",
                "yamet.{S}.{N}.{f}.det.out"),
    output:
        op.join(SIM_OUTPUT, "within_cell",
                "scmet.{S}.{N}.{f}.rds"),
    log:
        op.join("logs", "sim_within_cell_scmet_{S}_{N}_{f}.log"),
    threads: 4
    benchmark:
        repeat(op.join(SIM_OUTPUT, "within_cell",
                       "scmet.{S}.{N}.{f}.benchmark"), 5)
    script:
        op.join("src", "scmet.R")


rule sim_within_cell_plots:
    conda:
        op.join("..", "envs", "sim_r_plot.yaml")
    input:
        s_scmets=expand(
            op.join(SIM_OUTPUT, "within_cell",
                    "scmet.{S}.{N}.{f}.rds"),
            S=WC_S_GRID["S"], N=WC_S_GRID["N"], f=WC_S_GRID["f"],
        ),
        n_scmets=expand(
            op.join(SIM_OUTPUT, "within_cell",
                    "scmet.{S}.{N}.{f}.rds"),
            S=WC_N_GRID["S"], N=WC_N_GRID["N"], f=WC_N_GRID["f"],
        ),
        disp_scmet=op.join(
            SIM_OUTPUT, "within_cell",
            f"scmet.{WC_DISP['S']}.{WC_DISP['N']}.{WC_DISP['f']}.rds",
        ),
    params:
        S_GRID=WC_S_GRID,
        N_GRID=WC_N_GRID,
        DISP=WC_DISP,
    threads: 4
    output:
        html=op.join(SIM_RESULTS, "simulation_05_within_cell.html"),
    log:
        op.join("logs", "sim_within_cell_plots.log"),
    script:
        op.join("src", "simulations_05_within_cell.Rmd")


## -- between-cell variability simulations -------------------------------------

rule sim_between_cell_template:
    conda:
        op.join("..", "envs", "sim_r_sim.yaml")
    output:
        op.join(SIM_OUTPUT, "between_cell",
                "intervals.{S}.{N}.{f}.tsv"),
    log:
        op.join("logs", "sim_between_cell_template_{S}_{N}_{f}.log"),
    script:
        op.join("src", "simulations_between_cell_createTemplate.R")


rule sim_between_cell_samples:
    conda:
        op.join("..", "envs", "sim_r_sim.yaml")
    input:
        op.join(SIM_OUTPUT, "between_cell",
                "intervals.{S}.{N}.{f}.tsv"),
    output:
        op.join(SIM_OUTPUT, "between_cell",
                "sim.{s}.{S}.{N}.{f}.tsv"),
    log:
        op.join("logs", "sim_between_cell_samples_{s}_{S}_{N}_{f}.log"),
    params:
        data=op.join(SIM_OUTPUT, "between_cell"),
    script:
        op.join("src", "simulations_between_cell_getSamples.R")


rule sim_between_cell_yamet:
    conda:
        op.join("..", "envs", "sim_yamet.yaml")
    input:
        cells=lambda wildcards: expand(
            op.join(SIM_OUTPUT, "between_cell",
                    "sim.{s}.{{S}}.{{N}}.{{f}}.tsv"),
            s=range(1, int(wildcards.S) + 1),
        ),
        yamet=_SIM_YAMET,
    output:
        det_out=op.join(SIM_OUTPUT, "between_cell",
                        "yamet.{S}.{N}.{f}.det.out"),
        norm_det_out=op.join(SIM_OUTPUT, "between_cell",
                             "yamet.{S}.{N}.{f}.norm.det.out"),
    log:
        op.join("logs", "sim_between_cell_yamet_{S}_{N}_{f}.log"),
    params:
        data=op.join(SIM_OUTPUT, "between_cell"),
    threads: 4
    benchmark:
        repeat(op.join(SIM_OUTPUT, "between_cell",
                       "yamet.{S}.{N}.{f}.benchmark"), 5)
    shell:
        """
        {input.yamet} \
            --cell {params.data}/sim.*.{wildcards.S}.{wildcards.N}.{wildcards.f}.tsv \
            --reference {params.data}/ref.{wildcards.N}.{wildcards.f}.tsv \
            --intervals {params.data}/intervals.{wildcards.S}.{wildcards.N}.{wildcards.f}.tsv \
            --skip-header-intervals \
            --det-out {output.det_out} \
            --norm-det-out {output.norm_det_out} \
            --meth-out {params.data}/yamet.{wildcards.S}.{wildcards.N}.{wildcards.f}.meth.out \
            --cores {threads} \
            --no-print-sampens \
            &> {log}
        """


rule sim_between_cell_scmet:
    conda:
        op.join("..", "envs", "sim_r_plot.yaml")
    input:
        op.join(SIM_OUTPUT, "between_cell",
                "yamet.{S}.{N}.{f}.det.out"),
    output:
        op.join(SIM_OUTPUT, "between_cell",
                "scmet.{S}.{N}.{f}.rds"),
    log:
        op.join("logs", "sim_between_cell_scmet_{S}_{N}_{f}.log"),
    threads: 4
    benchmark:
        repeat(op.join(SIM_OUTPUT, "between_cell",
                       "scmet.{S}.{N}.{f}.benchmark"), 5)
    script:
        op.join("src", "scmet.R")


rule sim_between_cell_plots:
    conda:
        op.join("..", "envs", "sim_r_plot.yaml")
    input:
        s_scmets=expand(
            op.join(SIM_OUTPUT, "between_cell",
                    "scmet.{S}.{N}.{f}.rds"),
            S=BC_S_GRID["S"], N=BC_S_GRID["N"], f=BC_S_GRID["f"],
        ),
        n_scmets=expand(
            op.join(SIM_OUTPUT, "between_cell",
                    "scmet.{S}.{N}.{f}.rds"),
            S=BC_N_GRID["S"], N=BC_N_GRID["N"], f=BC_N_GRID["f"],
        ),
        disp_scmet=op.join(
            SIM_OUTPUT, "between_cell",
            f"scmet.{BC_DISP['S']}.{BC_DISP['N']}.{BC_DISP['f']}.rds",
        ),
    params:
        S_GRID=BC_S_GRID,
        N_GRID=BC_N_GRID,
        DISP=BC_DISP,
    threads: 4
    output:
        html=op.join(SIM_RESULTS, "simulation_06_between_cell.html"),
    log:
        op.join("logs", "sim_between_cell_plots.log"),
    script:
        op.join("src", "simulations_06_between_cell.Rmd")


## -- combined figures ---------------------------------------------------------

rule sim_combined_figure:
    conda:
        op.join("..", "envs", "sim_r_plot.yaml")
    input:
        within_html=op.join(SIM_RESULTS, "simulation_05_within_cell.html"),
        across_html=op.join(SIM_RESULTS, "simulation_06_between_cell.html"),
    output:
        html=op.join(SIM_RESULTS, "simulation_07_combined_figure.html"),
    log:
        op.join("logs", "sim_combined_figure.log"),
    script:
        op.join("src", "simulations_07_combined_figure.Rmd")


rule sim_combined_figure_adj:
    conda:
        op.join("..", "envs", "sim_r_plot.yaml")
    input:
        within_html=op.join(SIM_RESULTS, "simulation_05_within_cell.html"),
        across_html=op.join(SIM_RESULTS, "simulation_06_between_cell.html"),
    output:
        html=op.join(SIM_RESULTS, "simulation_08_combined_figure_adj.html"),
    log:
        op.join("logs", "sim_combined_figure_adj.log"),
    script:
        op.join("src", "simulations_08_combined_figure_adj.Rmd")
