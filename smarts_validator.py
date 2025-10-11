import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw, PandasTools
import io

# Page config
st.set_page_config(page_title="SMARTS Toolkit", layout="wide")

# Initialize session state
if 'mode' not in st.session_state:
    st.session_state.mode = "Visualizer"
if 'current_idx' not in st.session_state:
    st.session_state.current_idx = 0
if 'decisions' not in st.session_state:
    st.session_state.decisions = {}
if 'smarts_data' not in st.session_state:
    st.session_state.smarts_data = None
if 'test_molecules' not in st.session_state:
    st.session_state.test_molecules = None

# Title and mode selector
st.title("🔬 SMARTS Analysis Toolkit")

mode = st.radio(
    "Select Mode:",
    ["Visualizer", "Validator"],
    horizontal=True,
    help="Visualizer: Quick triage with traffic lights | Validator: Test against molecules"
)
st.session_state.mode = mode

st.write("---")

# ============================================================================
# MODE 1: VISUALIZER (Traffic Light Triage)
# ============================================================================

if mode == "Visualizer":
    st.subheader("📊 SMARTS Visualizer - Traffic Light Triage")
    
    uploaded_file = st.file_uploader(
        "Upload your SMARTS list (CSV with 'SMARTS' column)", 
        type=['csv'],
        key='visualizer_upload'
    )
    
    if uploaded_file:
        # Load SMARTS data
        if st.session_state.smarts_data is None:
            try:
                df = pd.read_csv(uploaded_file)
                if 'SMARTS' not in df.columns:
                    df.columns = ['SMARTS'] + list(df.columns[1:])
                st.session_state.smarts_data = df
                st.session_state.current_idx = 0
                st.session_state.decisions = {}
            except Exception as e:
                st.error(f"Could not parse file: {str(e)}")
                st.stop()
        
        df = st.session_state.smarts_data
        total = len(df)
        current = st.session_state.current_idx
        
        # Progress bar
        progress = len(st.session_state.decisions) / total if total > 0 else 0
        st.progress(progress)
        
        col_stat1, col_stat2, col_stat3, col_stat4, col_stat5 = st.columns(5)
        with col_stat1:
            st.metric("Progress", f"{len(st.session_state.decisions)}/{total}")
        with col_stat2:
            ok_count = sum(1 for d in st.session_state.decisions.values() if d == "OK")
            st.metric("OK", ok_count, delta_color="off")
        with col_stat3:
            yellow_count = sum(1 for d in st.session_state.decisions.values() if d == "Yellow")
            st.metric("Yellow", yellow_count, delta_color="off")
        with col_stat4:
            amber_count = sum(1 for d in st.session_state.decisions.values() if d == "Amber")
            st.metric("Amber", amber_count, delta_color="off")
        with col_stat5:
            red_count = sum(1 for d in st.session_state.decisions.values() if d == "Red")
            st.metric("Red", red_count, delta_color="off")
        
        st.write("---")
        
        # Display current SMARTS
        if current < total:
            col1, col2 = st.columns([2, 1])
            
            with col1:
                st.subheader(f"Pattern {current + 1} of {total}")
                smarts_pattern = df.iloc[current]['SMARTS']
                st.code(smarts_pattern, language='text')
                
                if 'Description' in df.columns and pd.notna(df.iloc[current]['Description']):
                    st.write(f"**Description:** {df.iloc[current]['Description']}")
                
                # Try to parse and display SMARTS
                try:
                    pattern = Chem.MolFromSmarts(smarts_pattern)
                    if pattern is None:
                        st.error("⚠️ INVALID SMARTS - Cannot parse!")
                    else:
                        st.success("✅ Valid SMARTS syntax")
                        
                        # Draw the SMARTS pattern itself
                        st.write("**Visual representation:**")
                        try:
                            img = Draw.MolToImage(pattern, size=(450, 350))
                            st.image(img, use_container_width=False)
                        except Exception as draw_error:
                            st.warning(f"Could not generate image: {str(draw_error)}")
                            
                except Exception as e:
                    st.error(f"Error processing SMARTS: {str(e)}")
            
            with col2:
                st.write("### Traffic Light Decision")
                
                # Show previous decision if exists
                if current in st.session_state.decisions:
                    prev_decision = st.session_state.decisions[current]
                    color_map = {
                        "OK": "🟢",
                        "Yellow": "🟡",
                        "Amber": "🟠",
                        "Red": "🔴"
                    }
                    st.info(f"Current: {color_map.get(prev_decision, '')} **{prev_decision}**")
                
                # Decision buttons (4 options)
                if st.button("🟢 OK", use_container_width=True, type="primary"):
                    st.session_state.decisions[current] = "OK"
                    if current < total - 1:
                        st.session_state.current_idx += 1
                    st.rerun()
                
                if st.button("🟡 Yellow", use_container_width=True):
                    st.session_state.decisions[current] = "Yellow"
                    if current < total - 1:
                        st.session_state.current_idx += 1
                    st.rerun()
                
                if st.button("🟠 Amber", use_container_width=True):
                    st.session_state.decisions[current] = "Amber"
                    if current < total - 1:
                        st.session_state.current_idx += 1
                    st.rerun()
                
                if st.button("🔴 Red", use_container_width=True):
                    st.session_state.decisions[current] = "Red"
                    if current < total - 1:
                        st.session_state.current_idx += 1
                    st.rerun()
                
                st.write("---")
                
                # Navigation
                col_prev, col_next = st.columns(2)
                with col_prev:
                    if st.button("⬅️ Previous", disabled=(current == 0)):
                        st.session_state.current_idx -= 1
                        st.rerun()
                with col_next:
                    if st.button("Next ➡️", disabled=(current >= total - 1)):
                        st.session_state.current_idx += 1
                        st.rerun()
                
                # Jump to specific index
                jump_to = st.number_input("Jump to #", min_value=1, max_value=total, 
                                          value=current + 1, step=1, key='jump_visualizer')
                if st.button("Go", key='go_visualizer'):
                    st.session_state.current_idx = jump_to - 1
                    st.rerun()
        
        else:
            st.success("🎉 You've reviewed all SMARTS patterns!")
        
        # Export results
        if len(st.session_state.decisions) > 0:
            st.write("---")
            st.subheader("📥 Export Results")
            
            # Create results dataframe
            results_df = df.copy()
            results_df['Decision'] = results_df.index.map(
                lambda x: st.session_state.decisions.get(x, "NOT_REVIEWED")
            )
            
            # Summary
            col1, col2, col3 = st.columns(3)
            with col1:
                reviewed = len([d for d in st.session_state.decisions.values() if d != "NOT_REVIEWED"])
                st.metric("Reviewed", reviewed)
            with col2:
                concerning = sum(1 for d in st.session_state.decisions.values() if d in ["Amber", "Red"])
                st.metric("Concerning (Amber/Red)", concerning)
            with col3:
                st.metric("Pending", total - len(st.session_state.decisions))
            
            # Download full results
            st.download_button(
                "📥 Download All Results",
                results_df.to_csv(index=False),
                "smarts_triage_results.csv",
                "text/csv",
                key='download_all'
            )
            
            # Download filtered lists
            col_d1, col_d2 = st.columns(2)
            with col_d1:
                concerning_df = results_df[results_df['Decision'].isin(['Amber', 'Red'])]
                if len(concerning_df) > 0:
                    st.download_button(
                        "📥 Download Concerning (Amber/Red)",
                        concerning_df.to_csv(index=False),
                        "concerning_smarts.csv",
                        "text/csv",
                        key='download_concerning'
                    )
            
            with col_d2:
                ok_df = results_df[results_df['Decision'] == 'OK']
                if len(ok_df) > 0:
                    st.download_button(
                        "📥 Download OK Patterns",
                        ok_df.to_csv(index=False),
                        "ok_smarts.csv",
                        "text/csv",
                        key='download_ok'
                    )
    
    else:
        st.info("👆 Upload a CSV file with your SMARTS patterns to begin triage")
        st.write("**Expected format:**")
        st.code("""SMARTS,Description
[#6]-[#8],Carbon-oxygen single bond
c1ccccc1,Benzene ring
[NX3;H2,H1;!$(NC=O)],Primary or secondary amine""")

# ============================================================================
# MODE 2: VALIDATOR (Test Against Molecules)
# ============================================================================

elif mode == "Validator":
    st.subheader("🧪 SMARTS Validator - Test Against Molecules")
    
    col_up1, col_up2 = st.columns(2)
    
    with col_up1:
        smarts_file = st.file_uploader(
            "Upload SMARTS patterns (CSV or SDF)",
            type=['csv', 'sdf'],
            key='validator_smarts'
        )
    
    with col_up2:
        molecules_file = st.file_uploader(
            "Upload test molecules (CSV with SMILES or SDF)",
            type=['csv', 'sdf'],
            key='validator_molecules'
        )
    
    # Load SMARTS
    if smarts_file and st.session_state.smarts_data is None:
        try:
            if smarts_file.name.endswith('.csv'):
                df = pd.read_csv(smarts_file)
                if 'SMARTS' not in df.columns:
                    df.columns = ['SMARTS'] + list(df.columns[1:])
            else:  # SDF
                df = PandasTools.LoadSDF(smarts_file)
                st.warning("SDF loaded - assuming 'ROMol' column contains patterns")
            st.session_state.smarts_data = df
            st.session_state.current_idx = 0
            st.session_state.decisions = {}
        except Exception as e:
            st.error(f"Error loading SMARTS: {str(e)}")
    
    # Load test molecules
    if molecules_file and st.session_state.test_molecules is None:
        try:
            if molecules_file.name.endswith('.csv'):
                mol_df = pd.read_csv(molecules_file)
                # Assume first column is SMILES or has 'SMILES' column
                if 'SMILES' in mol_df.columns:
                    smiles_col = 'SMILES'
                else:
                    smiles_col = mol_df.columns[0]
                
                mol_df['Mol'] = mol_df[smiles_col].apply(lambda x: Chem.MolFromSmiles(x) if pd.notna(x) else None)
                mol_df = mol_df[mol_df['Mol'].notna()]
                st.session_state.test_molecules = mol_df
                st.success(f"✅ Loaded {len(mol_df)} test molecules")
            else:  # SDF
                mol_df = PandasTools.LoadSDF(molecules_file)
                st.session_state.test_molecules = mol_df
                st.success(f"✅ Loaded {len(mol_df)} test molecules")
        except Exception as e:
            st.error(f"Error loading molecules: {str(e)}")
    
    # Main validator interface
    if st.session_state.smarts_data is not None and st.session_state.test_molecules is not None:
        df = st.session_state.smarts_data
        mol_df = st.session_state.test_molecules
        total = len(df)
        current = st.session_state.current_idx
        
        # Progress
        st.progress(current / total if total > 0 else 0)
        flagged_count = sum(1 for d in st.session_state.decisions.values() if d == "FLAGGED")
        col_m1, col_m2 = st.columns(2)
        with col_m1:
            st.metric("Current", f"{current + 1}/{total}")
        with col_m2:
            st.metric("Flagged", flagged_count, delta_color="off")
        
        st.write("---")
        
        if current < total:
            smarts_pattern = df.iloc[current]['SMARTS']
            
            col1, col2 = st.columns([2, 1])
            
            with col1:
                st.subheader(f"Testing Pattern {current + 1} of {total}")
                st.code(smarts_pattern, language='text')
                
                if 'Description' in df.columns and pd.notna(df.iloc[current]['Description']):
                    st.write(f"**Description:** {df.iloc[current]['Description']}")
                
                # Test against molecules
                try:
                    pattern = Chem.MolFromSmarts(smarts_pattern)
                    if pattern is None:
                        st.error("⚠️ INVALID SMARTS")
                    else:
                        # Find matches
                        matches = []
                        for idx, row in mol_df.iterrows():
                            mol = row['Mol'] if 'Mol' in row else row.get('ROMol')
                            if mol and mol.HasSubstructMatch(pattern):
                                matches.append((idx, mol))
                        
                        st.write(f"**Matches: {len(matches)} / {len(mol_df)} molecules ({len(matches)/len(mol_df)*100:.1f}%)**")
                        
                        if len(matches) > 0:
                            st.write("**Example matches (showing up to 6):**")
                            cols = st.columns(3)
                            for i, (idx, mol) in enumerate(matches[:6]):
                                with cols[i % 3]:
                                    try:
                                        img = Draw.MolToImage(mol, size=(200, 200), 
                                                             highlightAtoms=mol.GetSubstructMatch(pattern))
                                        st.image(img)
                                        st.caption(f"Mol {idx}")
                                    except:
                                        st.write(f"Mol {idx}")
                        else:
                            st.info("No matches found in test set")
                            
                except Exception as e:
                    st.error(f"Error: {str(e)}")
            
            with col2:
                st.write("### Actions")
                
                # Show if flagged
                if current in st.session_state.decisions:
                    if st.session_state.decisions[current] == "FLAGGED":
                        st.warning("🚩 **FLAGGED**")
                    else:
                        st.info("✓ Checked")
                
                # Flag button
                if st.button("🚩 FLAG THIS PATTERN", use_container_width=True, type="primary"):
                    st.session_state.decisions[current] = "FLAGGED"
                    st.rerun()
                
                if st.button("✓ Mark as Checked", use_container_width=True):
                    st.session_state.decisions[current] = "CHECKED"
                    if current < total - 1:
                        st.session_state.current_idx += 1
                    st.rerun()
                
                st.write("---")
                
                # Navigation
                col_prev, col_next = st.columns(2)
                with col_prev:
                    if st.button("⬅️ Previous", disabled=(current == 0), key='prev_val'):
                        st.session_state.current_idx -= 1
                        st.rerun()
                with col_next:
                    if st.button("Next ➡️", disabled=(current >= total - 1), key='next_val'):
                        st.session_state.current_idx += 1
                        st.rerun()
        
        else:
            st.success("🎉 Finished reviewing all patterns!")
        
        # Export
        if len(st.session_state.decisions) > 0:
            st.write("---")
            st.subheader("📥 Export Results")
            
            results_df = df.copy()
            results_df['Status'] = results_df.index.map(
                lambda x: st.session_state.decisions.get(x, "NOT_REVIEWED")
            )
            
            flagged_df = results_df[results_df['Status'] == 'FLAGGED']
            
            col_e1, col_e2 = st.columns(2)
            with col_e1:
                st.download_button(
                    "📥 Download All Results",
                    results_df.to_csv(index=False),
                    "validation_results.csv",
                    "text/csv",
                    key='val_download_all'
                )
            with col_e2:
                if len(flagged_df) > 0:
                    st.download_button(
                        "📥 Download Flagged Patterns Only",
                        flagged_df.to_csv(index=False),
                        "flagged_smarts.csv",
                        "text/csv",
                        key='val_download_flagged'
                    )
    
    else:
        st.info("👆 Upload both SMARTS patterns and test molecules to begin validation")
        st.write("**Note:** For CSV files, SMARTS should be in 'SMARTS' column and molecules in 'SMILES' column")
