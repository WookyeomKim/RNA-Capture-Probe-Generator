import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional
import matplotlib.pyplot as plt
from dataclasses import dataclass
import warnings
import os
import gc
import subprocess
from pathlib import Path
import time
import random

warnings.filterwarnings("ignore")

try:
    import nupack
    NUPACK_AVAILABLE = True
except ImportError:
    print("Warning: NUPACK not installed. Install with: pip install -U nupack")
    NUPACK_AVAILABLE = False


# ==========================================
# Global NUPACK Worker Function
# ==========================================
def nupack_worker_global(target_seq: str, probe_seq: str, temp_c: float, sodium_conc: float, a_conc: float, b_conc: float) -> float:
    if len(target_seq) < 10 or len(probe_seq) < 10:
        return 0.0

    try:
        model = nupack.Model(material="dna", celsius=temp_c, sodium=sodium_conc)
        strand_a = nupack.Strand(target_seq, name="a")
        strand_b = nupack.Strand(probe_seq, name="b")

        strands = {strand_a: a_conc, strand_b: b_conc}
        complex_ab = nupack.Complex([strand_a, strand_b])
        complexes = nupack.SetSpec(max_size=0, include=[complex_ab])

        my_set = nupack.ComplexSet(strands=strands, complexes=complexes)
        result = nupack.complex_analysis(my_set, compute=["pairs"], model=model)

        for key in result:
            energy = result[key].free_energy
            if np.isnan(energy) or np.isinf(energy) or energy > 0 or energy < -200:
                return 0.0
            return float(energy)

    except Exception:
        return 0.0
    return 0.0


@dataclass
class ProbeParameters:
    """
    ====================================================================
    💡 USER SETTINGS: Customize these parameters for your experiment
    ====================================================================
    """
    # ---------------------------------------------------------
    # [1] Input & Output Settings
    # ---------------------------------------------------------
    sequence_file: str = "sample.fasta"
    output_file_prefix: str = "capture_probes_ranked"
    max_sequences_to_process: Optional[int] = None
    
    # ---------------------------------------------------------
    # [2] Basic Probe Criteria
    # ---------------------------------------------------------
    probe_length: int = 70
    mismatch_tolerance: int = 9   
    min_gc_content: float = 40.0
    max_gc_content: float = 65.0
    
    # ---------------------------------------------------------
    # [3] NUPACK Thermodynamics Settings
    # ---------------------------------------------------------
    skip_nupack: bool = False
    deltaG_standard: float = -41  
    temp: float = 65.0            
    sodium_conc: float = 0.5      
    A_conc: float = 1.7e-17       
    B_conc: float = 4e-8          
    max_mismatch_for_nupack: int = 15 
    
    # ---------------------------------------------------------
    # [4] Accessibility Settings (RNAplfold)
    # ---------------------------------------------------------
    vienna_temp_c: float = 65.0   
    flank_nt: int = 150           
    punpaired_segment_len: int = 10 
    accessibility_threshold: float = 0 
    plot_accessibility: bool = True
    
    # ---------------------------------------------------------
    # [5] Computational Efficiency & Candidates Settings
    # ---------------------------------------------------------
    num_candidates: Optional[int] = None           
    candidate_range: float = 0.9                   
    batch_size: int = 50
    max_candidates_hard_limit: int = 5000          
    save_history_incrementally: bool = True
    history_chunk_size: int = 100
    temp_history_dir: str = "temp_history"
    save_full_deltaG: bool = False
    
    # ---------------------------------------------------------
    # [6] Visualization & Log Settings
    # ---------------------------------------------------------
    create_plots: bool = True
    save_plots: bool = False
    plot_dpi: int = 300
    verbose: bool = True
    log_invalid_segments: bool = True
    max_invalid_logs: int = 200          
    invalid_log_every: int = 1   


class MemoryOptimizedProbeDesigner:
    def __init__(self, params: ProbeParameters = None):
        self.params = params or ProbeParameters()
        self.iupac_map = self._create_iupac_map()
        self.results_history = []
        self.plot_counter = 0
        self.ini_match_results = None
        self.sequences_store = None
        self.consensus_store = None
        self.total_gc = None

        self.history_file_counter = 0
        self.temp_history_files = []
        self.structured_history = []

        self.access_array = None
        self.access_percentile = None
        self.low_access_mask = None
        
        self.global_nupack_cache = {}

        if self.params.save_history_incrementally:
            os.makedirs(self.params.temp_history_dir, exist_ok=True)

    def _create_iupac_map(self) -> Dict[str, str]:
        return {
            "A": "A", "C": "C", "G": "G", "T": "T",
            "R": "AG", "Y": "CT", "S": "CG", "W": "AT",
            "K": "GT", "M": "AC", "B": "CGT", "D": "AGT",
            "H": "ACT", "V": "ACG", "N": "ACGT", "-": "-"
        }

    def save_history_chunk(self):
        if not self.structured_history:
            return
        chunk_file = os.path.join(self.params.temp_history_dir, f"history_chunk_{self.history_file_counter:04d}.pkl")
        import pickle
        with open(chunk_file, "wb") as f:
            pickle.dump(self.structured_history, f)
        self.temp_history_files.append(chunk_file)
        self.history_file_counter += 1
        self.structured_history = []
        gc.collect()

    def load_all_history_chunks(self):
        import pickle
        full_history = []
        for chunk_file in self.temp_history_files:
            with open(chunk_file, "rb") as f:
                chunk = pickle.load(f)
                full_history.extend(chunk)
        full_history.extend(self.structured_history)
        return full_history

    def cleanup_temp_files(self):
        import shutil
        if os.path.exists(self.params.temp_history_dir):
            shutil.rmtree(self.params.temp_history_dir)

    def read_fasta(self, filename: str) -> List[str]:
        sequences = []
        current_seq = []
        with open(filename, "r") as f:
            for line in f:
                line = line.strip()
                if not line: continue
                if line.startswith(">"):
                    if current_seq:
                        sequences.append("".join(current_seq))
                        current_seq = []
                else:
                    current_seq.append(line)
            if current_seq:
                sequences.append("".join(current_seq))
        return sequences

    def generate_consensus(self, sequences: List[str], gap_threshold: float = 0.5) -> str:
        if not sequences: return ""
        alignment_length = len(sequences[0])
        num_seqs = len(sequences)
        nucleotides = "ACGT"
        nucleotide_map = {"A": 0, "C": 1, "G": 2, "T": 3}

        counts = np.zeros((4, alignment_length), dtype=int)
        gap_counts = np.zeros(alignment_length, dtype=int)

        for seq in sequences:
            for j, base in enumerate(seq[:alignment_length]):
                if base in nucleotide_map:
                    counts[nucleotide_map[base], j] += 1
                elif base == "-":
                    gap_counts[j] += 1

        gap_ratio = gap_counts / num_seqs
        max_indices = np.argmax(counts, axis=0)
        total_acgt = np.sum(counts, axis=0)

        consensus = []
        for j in range(alignment_length):
            if gap_ratio[j] >= gap_threshold or total_acgt[j] == 0:
                consensus.append("-")
            else:
                consensus.append(nucleotides[max_indices[j]])
        return "".join(consensus)

    def calculate_gc_content(self, sequence: str) -> float:
        seq_clean = sequence.replace("-", "")
        if not seq_clean: return 0.0
        gc_count = seq_clean.count("G") + seq_clean.count("C")
        return (gc_count / len(seq_clean)) * 100.0

    def reverse_complement_iupac(self, sequence: str) -> str:
        complement_map = {
            "A": "T", "T": "A", "G": "C", "C": "G",
            "R": "Y", "Y": "R", "S": "S", "W": "W",
            "K": "M", "M": "K", "B": "V", "V": "B",
            "D": "H", "H": "D", "N": "N", "-": "-"
        }
        seq_clean = sequence.replace("-", "")
        return "".join(complement_map.get(base, base) for base in seq_clean[::-1])

    def randomize_iupac(self, sequence: str) -> str:
        iupac_options = {
            "A": ["A"], "T": ["T"], "G": ["G"], "C": ["C"],
            "R": ["A", "G"], "Y": ["C", "T"], "S": ["G", "C"], "W": ["A", "T"],
            "K": ["G", "T"], "M": ["A", "C"], "B": ["C", "G", "T"], "D": ["A", "G", "T"],
            "H": ["A", "C", "T"], "V": ["A", "C", "G"], "N": ["A", "T", "G", "C"],
            "-": []
        }
        result = []
        for base in sequence:
            if base in iupac_options:
                options = iupac_options[base]
                if options: result.append(random.choice(options))
            elif base in "ATGC":
                result.append(base)
            else:
                result.append("A")
        return "".join(result)

    def create_match_matrix(self, sequences: List[str], consensus: str) -> np.ndarray:
        num_seqs = len(sequences)
        seq_length = len(consensus)
        match_matrix = np.zeros((num_seqs, seq_length), dtype=np.int8)

        for i, seq in enumerate(sequences):
            if i % 1000 == 0 and i > 0 and self.params.verbose:
                print(f"Constructing match matrix: {i}")

            for j in range(min(len(seq), seq_length)):
                if consensus[j] in self.iupac_map:
                    options = self.iupac_map[consensus[j]]
                    match_matrix[i, j] = 1 if (j < len(seq) and seq[j] in options) else 0
                else:
                    match_matrix[i, j] = 1 if (j < len(seq) and seq[j] == consensus[j]) else 0
        return match_matrix

    def find_gaps_in_probe_regions(self, consensus: str, probe_length: int) -> np.ndarray:
        seq_length = len(consensus)
        num_gaps = np.full(seq_length, -1, dtype=int)

        for i in range(seq_length):
            if i >= seq_length - probe_length + 1: break
            if consensus[i] == "-": continue

            actual_bases = 0
            j = i
            while j < seq_length and actual_bases < probe_length:
                if consensus[j] != "-": actual_bases += 1
                j += 1

            if actual_bases == probe_length:
                num_gaps[i] = j - i - probe_length
        return num_gaps

    def calculate_mismatch_tolerance_with_record(self, match_matrix: np.ndarray, consensus: str, num_gaps: np.ndarray, mismatch_tolerance: int):
        seq_length = len(consensus)
        num_sequences = match_matrix.shape[0]
        probe_length = self.params.probe_length

        valid_positions = seq_length - probe_length + 1
        tolerance_results = np.full(valid_positions, -np.inf)
        gc_contents = np.full(valid_positions, np.nan)

        i = 0
        while i < valid_positions:
            if i < len(consensus) and consensus[i] == "-":
                i += 1
            else:
                if i < len(num_gaps) and num_gaps[i] != -1:
                    end_index = min(i + probe_length + num_gaps[i], seq_length)
                    if end_index <= seq_length:
                        probe_matches = match_matrix[:, i:end_index]
                        mismatch_counts = np.sum(probe_matches == 0, axis=1)
                        tolerance_results[i] = np.sum(mismatch_counts <= mismatch_tolerance)
                        
                        current_seq = consensus[i:end_index]
                        gc_contents[i] = self.calculate_gc_content(current_seq)
                i += 1
        return tolerance_results, gc_contents, None

    def visualize_results(self, raw_tolerance, gc_contents, candidate_positions=None, save_plots=False):
        min_length = min(len(raw_tolerance), len(gc_contents))
        raw_tolerance = raw_tolerance[:min_length]
        gc_contents = gc_contents[:min_length]

        print("\n📊 [1/4] Figure 1: Mismatch Tolerance Distribution (Before filtering)")
        plt.figure(figsize=(10, 6))
        plt.plot(raw_tolerance)
        plt.xlabel("Position (bp)", fontsize=14)
        plt.ylabel("Number of capturable targets", fontsize=14)
        plt.title("Criteria: Mismatch tolerance (Before filtering)", fontsize=14)
        plt.grid(True, alpha=0.3)
        if save_plots: plt.savefig(f"figure1_raw_tolerance_{self.plot_counter}.png", dpi=self.params.plot_dpi)
        plt.show()

        print("📊 [2/4] Figure 2: GC Content Distribution")
        plt.figure(figsize=(10, 6))
        plt.plot(gc_contents)
        plt.xlabel("Position (bp)", fontsize=14)
        plt.ylabel("GC contents(%)", fontsize=14)
        plt.title("GC contents over the whole sequence", fontsize=14)
        plt.grid(True, alpha=0.3)
        if save_plots: plt.savefig(f"figure2_gc_content_{self.plot_counter}.png", dpi=self.params.plot_dpi)
        plt.show()

        filtered_tolerance = raw_tolerance.copy()
        
        limit_gc = (gc_contents >= self.params.min_gc_content) & (gc_contents <= self.params.max_gc_content)
        filtered_tolerance[~limit_gc] = np.nan
        
        if self.low_access_mask is not None:
            n = min(len(filtered_tolerance), len(self.low_access_mask))
            filtered_tolerance[:n][self.low_access_mask[:n]] = np.nan

        print("📊 [3/4] Figure 3: Tolerance after GC & Accessibility Filtering")
        plt.figure(figsize=(10, 6))
        plt.plot(filtered_tolerance)
        plt.xlabel("Position (bp)", fontsize=14)
        plt.ylabel("Number of capturable targets", fontsize=14)
        plt.title("Tolerance after GC & Accessibility Filtering", fontsize=14)
        plt.grid(True, alpha=0.3)
        if save_plots: plt.savefig(f"figure3_filtered_tolerance_{self.plot_counter}.png", dpi=self.params.plot_dpi)
        plt.show()

        if candidate_positions is not None and len(candidate_positions) > 0:
            print(f"📊 [4/4] Figure 4: Final candidate positions (n={len(candidate_positions)})")
            plt.figure(figsize=(10, 6))
            plt.plot(filtered_tolerance, alpha=0.7)
            candidate_values = filtered_tolerance[candidate_positions]
            plt.scatter(candidate_positions, candidate_values, color="red", s=30, zorder=5)
            plt.xlabel("Position (bp)", fontsize=14)
            plt.ylabel("Number of capturable targets", fontsize=14)
            plt.title("Criteria: Mismatch tolerance + GC contents (with candidates)", fontsize=14)
            plt.grid(True, alpha=0.3)
            if save_plots: plt.savefig(f"figure4_candidates_{self.plot_counter}.png", dpi=self.params.plot_dpi)
            plt.show()

        self.plot_counter += 1

    def plot_accessibility(self, access_array: np.ndarray, low_mask: np.ndarray, threshold: float):
        print("\n📊 Accessibility Scan Results")
        x = np.arange(1, len(access_array) + 1)
        y = access_array.copy()

        plt.figure(figsize=(12, 5))
        plt.plot(x, y, alpha=0.8)
        plt.axhline(threshold, color="red", linestyle=":", linewidth=1)

        low_x = x[low_mask & ~np.isnan(y)]
        low_y = y[low_mask & ~np.isnan(y)]
        if len(low_x) > 0:
            plt.scatter(low_x, low_y, color="red", s=12, zorder=5)

        plt.xlabel("Probe start position in gapped consensus (1-based)", fontsize=12)
        plt.ylabel("Accessibility mean P(unpaired len=10) in 70 nt window", fontsize=12)
        plt.title("Accessibility scan for all 70 nt windows", fontsize=12)
        plt.grid(True, alpha=0.3)
        plt.show()

    def _run_rnaplfold_u10_local(self, rna_seq: str, prefix: str, W: int, L: int) -> Path:
        rna_seq = rna_seq.strip().upper().replace("T", "U").replace("-", "")
        fasta = f">{prefix}\n{rna_seq}\n"
        cmd = ["RNAplfold", "-u", str(self.params.punpaired_segment_len), "-W", str(W), "-L", str(L), "-T", str(self.params.vienna_temp_c)]
        subprocess.run(cmd, input=fasta.encode(), check=True)
        return Path(f"{prefix}_lunp")

    def _load_p10_start_from_lunp(self, lunp_path: str) -> pd.Series:
        seg_len = self.params.punpaired_segment_len
        col_index = seg_len
        starts, vals = [], []
        with open(lunp_path, "r", encoding="utf-8") as f:
            for line in f:
                s = line.strip()
                if (not s) or s.startswith("#"): continue
                parts = s.split()
                if len(parts) <= col_index: continue
                i_end = int(parts[0])
                v = parts[col_index]
                if v == "NA": continue
                i_start = i_end - seg_len + 1
                if i_start < 1: continue
                starts.append(i_start)
                vals.append(float(v))
        return pd.Series(vals, index=pd.Index(starts, name="i_start"), name=f"p{seg_len}_start").sort_index()

    def compute_accessibility_scan_and_mask(self, consensus_gapped: str) -> Tuple[np.ndarray, np.ndarray]:
        probe_len = self.params.probe_length
        flank = self.params.flank_nt
        seg_len = self.params.punpaired_segment_len
        threshold = self.params.accessibility_threshold

        seq_len_gapped = len(consensus_gapped)
        valid_positions = seq_len_gapped - probe_len + 1

        nongap_prefix = np.zeros(seq_len_gapped + 1, dtype=int)
        for i, ch in enumerate(consensus_gapped):
            nongap_prefix[i + 1] = nongap_prefix[i] + (1 if ch != "-" else 0)

        ungapped = consensus_gapped.replace("-", "")
        W = probe_len + 2 * flank
        L = W
        prefix = "consensus_access"
        lunp = self._run_rnaplfold_u10_local(ungapped, prefix=prefix, W=W, L=L)
        p_start = self._load_p10_start_from_lunp(str(lunp))

        num_gaps = self.find_gaps_in_probe_regions(consensus_gapped, probe_len)
        access_array = np.full(valid_positions, np.nan, dtype=float)
        low_mask = np.zeros(valid_positions, dtype=bool)

        for i in range(valid_positions):
            if consensus_gapped[i] == "-": continue
            if i >= len(num_gaps) or num_gaps[i] == -1: continue
            
            end_index = i + probe_len + num_gaps[i]
            if end_index > seq_len_gapped: continue

            start_ungapped_1based = nongap_prefix[i] + 1
            last_start = start_ungapped_1based + (probe_len - seg_len)
            window_vals = p_start.loc[start_ungapped_1based:last_start]
            if len(window_vals) == 0: continue

            a = float(window_vals.mean())
            access_array[i] = a
            if a < threshold: low_mask[i] = True

        if self.params.plot_accessibility:
            self.plot_accessibility(access_array, low_mask, threshold)

        valid = ~np.isnan(access_array)
        access_percentile = np.full_like(access_array, np.nan, dtype=float)
        if np.any(valid):
            s = pd.Series(access_array[valid])
            access_percentile[valid] = s.rank(pct=True, method="average").values * 100.0
        self.access_percentile = access_percentile

        return access_array, low_mask

    def find_optimal_probe_memory_optimized(self, sequences: List[str], consensus: str, mismatch_tolerance: int) -> None:
        num_sequences = len(sequences)
        seq_length = len(consensus)

        match_matrix = self.create_match_matrix(sequences, consensus)
        num_gaps = self.find_gaps_in_probe_regions(consensus, self.params.probe_length)
        
        tolerance_results, gc_contents, _ = self.calculate_mismatch_tolerance_with_record(
            match_matrix, consensus, num_gaps, mismatch_tolerance
        )
        gc.collect()

        print(f"\n🔍 Filtering Summary (MM <= {mismatch_tolerance})")
        tolerance_results_filtered = tolerance_results.copy()
        
        limit_gc = (gc_contents >= self.params.min_gc_content) & (gc_contents <= self.params.max_gc_content)
        gc_fail_count = np.sum(~limit_gc)
        print(f" ❌ GC Filter ({self.params.min_gc_content}~{self.params.max_gc_content}%): {gc_fail_count} positions dropped")
        tolerance_results_filtered[~limit_gc] = -np.inf
        tolerance_results_filtered[np.isnan(tolerance_results_filtered)] = -np.inf

        if self.low_access_mask is not None:
            min_len = min(len(tolerance_results_filtered), len(self.low_access_mask))
            access_fail_count = np.sum(self.low_access_mask[:min_len])
            print(f" ❌ Accessibility Filter (<{self.params.accessibility_threshold}): {access_fail_count} positions dropped")
            tolerance_results_filtered[:min_len][self.low_access_mask[:min_len]] = -np.inf

        sorted_indices = np.argsort(tolerance_results_filtered)[::-1]
        sorted_tolerance = tolerance_results_filtered[sorted_indices]

        max_tolerance = np.max(tolerance_results_filtered)
        if max_tolerance == -np.inf:
            print("No valid candidates passed the filters.")
            return

        candidate_cut = max_tolerance * self.params.candidate_range
        num_candidates = np.sum(tolerance_results_filtered >= candidate_cut)

        print(f" ✅ Final Candidates passing all filters: {num_candidates}")
        print("-" * 60)

        if num_candidates > self.params.max_candidates_hard_limit:
            num_candidates = self.params.max_candidates_hard_limit
        if self.params.num_candidates is not None and num_candidates > self.params.num_candidates:
            num_candidates = self.params.num_candidates

        top_values = sorted_tolerance[:num_candidates]
        candidate_positions = sorted_indices[:num_candidates]

        if self.params.create_plots:
            self.visualize_results(tolerance_results, gc_contents, candidate_positions, self.params.save_plots)

        batch_size = self.params.batch_size
        total_batches = (num_candidates - 1) // batch_size + 1

        loop_start_time = time.time()

        for batch_start in range(0, num_candidates, batch_size):
            batch_idx = batch_start // batch_size
            batch_end = min(batch_start + batch_size, num_candidates)
            batch_positions = candidate_positions[batch_start:batch_end]
            batch_top_values = top_values[batch_start:batch_end]

            if self.params.verbose:
                print("\n" + "=" * 60)
                print(f"▶️ [Batch {batch_idx + 1}/{total_batches}] Processing Candidates {batch_start + 1} ~ {batch_end}")

            batch_hyb_g = np.zeros((num_sequences, len(batch_positions)))
            batch_seq_history = []

            for m_batch, pos in enumerate(batch_positions):
                candidate_position_end = min(pos + self.params.probe_length + num_gaps[pos], seq_length)
                candidate_probe_sequence = consensus[pos:candidate_position_end]
                
                candidate_probe_clean = candidate_probe_sequence.replace("-", "")
                candidate_probe_clean_rc = self.reverse_complement_iupac(candidate_probe_clean)
                candidate_probe_fixed = self.randomize_iupac(candidate_probe_clean_rc)

                if not self.params.skip_nupack:
                    processed_targets = []
                    mismatch_exceed_count = 0
                    
                    for s, current_sequence in enumerate(sequences):
                        if pos < len(current_sequence):
                            candidate_target_sequence = current_sequence[pos:min(candidate_position_end, len(current_sequence))]
                            target_clean = candidate_target_sequence.replace("-", "")
                            
                            if len(target_clean) < 10:
                                processed_targets.append("TOO_SHORT")
                            else:
                                min_l = min(len(target_clean), len(candidate_probe_clean))
                                mm_count = sum(1 for a, b in zip(target_clean[:min_l], candidate_probe_clean[:min_l]) if a != b)
                                
                                if mm_count > self.params.max_mismatch_for_nupack:
                                    processed_targets.append("MISMATCH_EXCEED")
                                    mismatch_exceed_count += 1
                                else:
                                    target_final = self.randomize_iupac(target_clean)
                                    processed_targets.append(target_final)
                        else:
                            processed_targets.append("TOO_SHORT")
                    
                    unique_targets = list(set(processed_targets))
                    for trash in ["TOO_SHORT", "MISMATCH_EXCEED"]:
                        if trash in unique_targets:
                            unique_targets.remove(trash)
                    
                    if unique_targets:
                        total_u = len(unique_targets)
                        print(f"   [Candidate {batch_start + m_batch + 1}/{num_candidates}] Skipped: {mismatch_exceed_count} | Unique Targets: {total_u} ... Starting Cache Eval...")
                        
                        completed = 0
                        print_interval = max(1, total_u // 10)
                        
                        for seq in unique_targets:
                            cache_key = (seq, candidate_probe_fixed)
                            
                            if cache_key not in self.global_nupack_cache:
                                self.global_nupack_cache[cache_key] = nupack_worker_global(
                                    seq, candidate_probe_fixed, 
                                    self.params.temp, self.params.sodium_conc, 
                                    self.params.A_conc, self.params.B_conc
                                )
                            
                            completed += 1
                            if completed % print_interval == 0 or completed == total_u:
                                pct = (completed / total_u) * 100
                                print(f"      ... NUPACK Progress: {completed} / {total_u} ({pct:.1f}%)")

                    for s, tgt in enumerate(processed_targets):
                        if tgt in ["TOO_SHORT", "MISMATCH_EXCEED"]:
                            batch_hyb_g[s, m_batch] = 0.0
                        else:
                            batch_hyb_g[s, m_batch] = self.global_nupack_cache[(tgt, candidate_probe_fixed)]
                else:
                    batch_hyb_g[:, m_batch] = np.random.uniform(-110, -80, num_sequences)

                batch_seq_history.append(candidate_probe_fixed)

            batch_capture_result = (batch_hyb_g <= self.params.deltaG_standard).astype(int)
            batch_capture_counts = np.sum(batch_capture_result, axis=0)

            for i in range(len(batch_positions)):
                pos0 = int(batch_positions[i])
                candidate_data = {
                    "Initial_Matches_(Mismatch)": int(batch_top_values[i]),
                    "Position": pos0 + 1,
                    "Sequence": batch_seq_history[i],
                    "Captured_Targets_(DeltaG)": int(batch_capture_counts[i]),
                    "Coverage_Percent": round((float(batch_capture_counts[i]) / num_sequences) * 100.0, 2) if num_sequences > 0 else 0.0
                }
                
                if self.access_array is not None and pos0 < len(self.access_array):
                    candidate_data["Accessibility_Mean"] = round(float(self.access_array[pos0]), 4)
                else:
                    candidate_data["Accessibility_Mean"] = "N/A"
                    
                if self.access_percentile is not None and pos0 < len(self.access_percentile):
                    candidate_data["Accessibility_Percentile(%)"] = round(float(self.access_percentile[pos0]), 2)
                else:
                    candidate_data["Accessibility_Percentile(%)"] = "N/A"

                if not self.params.save_full_deltaG:
                    candidate_data.update({
                        "Avg_DeltaG": round(float(np.mean(batch_hyb_g[:, i])), 2), 
                        "Std_DeltaG": round(float(np.std(batch_hyb_g[:, i])), 2),
                        "Min_DeltaG": round(float(np.min(batch_hyb_g[:, i])), 2), 
                        "Max_DeltaG": round(float(np.max(batch_hyb_g[:, i])), 2)
                    })
                self.structured_history.append(candidate_data)

            if self.params.save_history_incrementally and len(self.structured_history) >= self.params.history_chunk_size:
                self.save_history_chunk()

            del batch_hyb_g
            del batch_capture_result
            gc.collect()

            if self.params.verbose:
                elapsed = time.time() - loop_start_time
                avg_time = elapsed / (batch_idx + 1)
                eta_secs = (total_batches - (batch_idx + 1)) * avg_time
                print(f" 🏁 [Batch {batch_idx + 1} completed] Elapsed time: {int(elapsed//60)}m {int(elapsed%60)}s | ETA: {int(eta_secs // 60)}m {int(eta_secs % 60)}s\n")

    def run_single_design(self, sequences: List[str], consensus: str) -> List[Dict]:
        self.sequences_store = sequences.copy()
        self.consensus_store = consensus

        self.find_optimal_probe_memory_optimized(
            self.sequences_store, self.consensus_store, self.params.mismatch_tolerance
        )
        
        full_history = self.load_all_history_chunks() if self.params.save_history_incrementally else self.structured_history
        
        # Sort by final NUPACK capture count to generate ranking
        sorted_history = sorted(full_history, key=lambda x: x.get('Captured_Targets_(DeltaG)', 0), reverse=True)
        
        # Add Rank
        for rank, item in enumerate(sorted_history):
            item["Rank"] = rank + 1
            
        # Reorder columns to match requested output format perfectly
        final_list = []
        for item in sorted_history:
            ordered_item = {
                "Rank": item.get("Rank"),
                "Initial_Matches_(Mismatch)": item.get("Initial_Matches_(Mismatch)"),
                "Captured_Targets_(DeltaG)": item.get("Captured_Targets_(DeltaG)"),
                "Coverage_Percent": item.get("Coverage_Percent"),
                "Position": item.get("Position"),
                "Sequence": item.get("Sequence"),
                "Accessibility_Mean": item.get("Accessibility_Mean"),
                "Accessibility_Percentile(%)": item.get("Accessibility_Percentile(%)"),
                "Avg_DeltaG": item.get("Avg_DeltaG"),
                "Std_DeltaG": item.get("Std_DeltaG"),
                "Min_DeltaG": item.get("Min_DeltaG"),
                "Max_DeltaG": item.get("Max_DeltaG")
            }
            final_list.append(ordered_item)

        return final_list

    def save_results_history_vertical(self, sorted_history: List[Dict], filename: str):
        if not sorted_history: 
            return pd.DataFrame()
        
        df_history = pd.DataFrame(sorted_history)
        df_history.to_excel(filename, index=False)
        
        if self.params.save_history_incrementally: 
            self.cleanup_temp_files()
        return df_history


def main(custom_params=None):
    params = custom_params or ProbeParameters()
    designer = MemoryOptimizedProbeDesigner(params)

    try:
        print(f"Reading sequence data from: {params.sequence_file}")
        sequences_all = designer.read_fasta(params.sequence_file)
        print("Sequences data has been received")
        print(f"Loaded {len(sequences_all)} sequences (full set)")

        consensus = designer.generate_consensus(sequences_all)
        print("Consensus generated successfully from input sequences.")

        if params.max_sequences_to_process is not None:
            sequences = sequences_all[:params.max_sequences_to_process]
            print(f"Limited to {params.max_sequences_to_process} sequences for testing")
        else:
            sequences = sequences_all

        print(f"Loaded {len(sequences)} sequences for probe design")
        
        print("\nStarting accessibility scan using RNAplfold (u=10) across all 70 nt windows")
        print(f"Flank length on each side: {params.flank_nt} nt")
        print(f"Accessibility threshold: {params.accessibility_threshold}")

        access_array, low_mask = designer.compute_accessibility_scan_and_mask(consensus)
        designer.access_array = access_array
        designer.low_access_mask = low_mask

        n_low = int(np.sum(low_mask))
        print(f"Accessibility scan done. Low accessibility windows: {n_low}")

        print("\n" + "=" * 50)
        print("🚀 Starting fast memory-optimized probe design (Single Pass)")
        print("=" * 50)

        ranked_probes = designer.run_single_design(sequences, consensus)
        
        if ranked_probes:
            output_name = f"{params.output_file_prefix}_MM{params.mismatch_tolerance}_ranked.xlsx"
            designer.save_results_history_vertical(ranked_probes, output_name)

            print("\n" + "=" * 50)
            print(f"✅ Data successfully saved to Excel: {output_name}")
            print(f"Total candidates evaluated: {len(ranked_probes)}")
            print(f"🏆 Best Probe (#1) covers {ranked_probes[0]['Coverage_Percent']}% of targets.")
            print("=" * 50)
        else:
            print("\n⚠️ No viable probes found matching the criteria.")
        
        return ranked_probes

    except Exception as e:
        print(f"Error occurred: {e}")
        import traceback
        traceback.print_exc()
        return None

    finally:
        if hasattr(designer, "params") and designer.params.save_history_incrementally:
            designer.cleanup_temp_files()

if __name__ == "__main__":
    main()