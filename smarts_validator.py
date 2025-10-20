import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw, PandasTools
import io
import requests
import base64
import time

# Page config
st.set_page_config(page_title="SMARTS Toolkit", layout="wide", initial_sidebar_state="collapsed")

# Custom CSS for compact layout
st.markdown("""
<style>
    .block-container {
        padding-top: 1rem;
        padding-bottom: 0rem;
        padding-left: 2rem;
        padding-right: 2rem;
    }
    [data-testid="stMetric"] {
        background-color: #f0f2f6;
        padding: 0.3rem 0.5rem;
        border-radius: 0.3rem;
    }
    [data-testid="stMetricLabel"] {
        font-size: 0.75rem;
    }
    [data-testid="stMetricValue"] {
        font-size: 1rem;
    }
    .stButton button {
        padding: 0.25rem 0.5rem;
        font-size: 0.85rem;
    }
    h1 {
        font-size: 1.5rem;
        margin-bottom: 0.5rem;
    }
    h2 {
        font-size: 1.2rem;
        margin-bottom: 0.5rem;
    }
    h3 {
        font-size: 1rem;
        margin-bottom: 0.3rem;
    }
    [data-testid="stFileUploader"] {
        padding: 0.5rem;
    }
    .element-container {
        margin-bottom: 0.5rem;
    }
</style>
""", unsafe_allow_html=True)

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
if 'api_key' not in st.session_state:
    st.session_state.api_key = ""
if 'viz_cache' not in st.session_state:
    st.session_state.viz_cache = {}
if 'filter_results' not in st.session_state:
    st.session_state.filter_results = None
if 'uploaded_file_name' not in st.session_state:
    st.session_state.uploaded_file_name = None

# Helper function for SMARTS.plus visualization
def get_smartsplus_image(smarts_pattern, api_key, use_cache=True):
    """Get visualization from SMARTS.plus API - returns SVG string"""
    cache_key = f"{smarts_pattern}_{api_key[:8]}_svg"
    if use_cache and cache_key in st.session_state.viz_cache:
        return st.session_state.viz_cache[cache_key]
    
    try:
        payload = {
            "query": {
                "smarts": smarts_pattern,
                "parameters": {
                    "file_format": "svg",  # Use SVG for vector quality
                    "visualization_mode": 0,
                    "legend_mode": 1,
                    "smarts_string_into_picture": True,
                    "visualization_of_default_bonds": 0
                }
            }
        }
        
        headers = {
            'Content-Type': 'application/json',
            'X-API-Key': api_key
        }
        
        response = requests.post(
            'https://api.smarts.plus/smartsView/',
            json=payload,
            headers=headers,
            timeout=10
        )
        
        if response.status_code == 202:
            job_id = response.json().get('job_id')
            for _ in range(5):
                time.sleep(1)
                result_response = requests.get(
                    f'https://api.smarts.plus/smartsView/?job_id={job_id}',
                    timeout=10
                )
                if result_response.status_code == 200:
                    result = result_response.json()
                    if 'result' in result and 'image' in result['result']:
                        svg_data = result['result']['image']
                        st.session_state.viz_cache[cache_key] = svg_data
                        return svg_data
        elif response.status_code == 200:
            result = response.json()
            if 'result' in result and 'image' in result['result']:
                svg_data = result['result']['image']
                st.session_state.viz_cache[cache_key] = svg_data
                return svg_data
        return None
    except Exception as e:
        # Show error in debug mode
        if st.session_state.get('debug_mode', False):
            st.error(f"SMARTS.plus error: {str(e)}")
        return None

# Header
st.title("🔬 SMARTS Toolkit")
mode = st.radio("Select Mode:", ["Visualizer", "Validator", "Gen AI Filter"], horizontal=True, key="mode_selector")
st.session_state.mode = mode

# Settings expander
with st.expander("⚙️ Settings", expanded=False):
    col_s1, col_s2 = st.columns(2)
    with col_s1:
        api_key_input = st.text_input(
            "SMARTS.plus API Key",
            type="password",
            value=st.session_state.api_key,
            help="Get free key at smarts.plus/sign_up"
        )
        if api_key_input != st.session_state.api_key:
            st.session_state.api_key = api_key_input
            st.session_state.viz_cache = {}
    with col_s2:
        use_smartsplus = st.checkbox(
            "Use SMARTS.plus viz",
            value=bool(st.session_state.api_key),
            disabled=not bool(st.session_state.api_key)
        )
        st.session_state.debug_mode = st.checkbox(
            "Debug mode",
            value=False,
            help="Show SMARTS.plus API errors"
        )

st.divider()

# ============================================================================
# MODE 1: VISUALIZER
# ============================================================================

if mode == "Visualizer":
    uploaded_file = st.file_uploader("📁 Upload SMARTS CSV", type=['csv'], label_visibility="collapsed")
    
    if uploaded_file:
        # Detect new file uploads
        current_file_id = f"{uploaded_file.name}_{uploaded_file.size}"
        
        if st.session_state.uploaded_file_name != current_file_id:
            try:
                df = pd.read_csv(uploaded_file)
                if 'SMARTS' not in df.columns:
                    df.columns = ['SMARTS'] + list(df.columns[1:])
                st.session_state.smarts_data = df
                st.session_state.current_idx = 0
                st.session_state.decisions = {}
                st.session_state.uploaded_file_name = current_file_id
            except Exception as e:
                st.error(f"Parse error: {str(e)}")
                st.stop()
        
        df = st.session_state.smarts_data
        total = len(df)
        current = st.session_state.current_idx
        progress = len(st.session_state.decisions) / total if total > 0 else 0
        
        st.progress(progress, text=f"Pattern {current + 1}/{total}")
        
        col1, col2, col3, col4, col5 = st.columns(5)
        with col1:
            st.metric("✓", f"{len(st.session_state.decisions)}/{total}", label_visibility="collapsed")
        with col2:
            st.metric("🟢", sum(1 for d in st.session_state.decisions.values() if d == "OK"), label_visibility="collapsed")
        with col3:
            st.metric("🟡", sum(1 for d in st.session_state.decisions.values() if d == "Yellow"), label_visibility="collapsed")
        with col4:
            st.metric("🟠", sum(1 for d in st.session_state.decisions.values() if d == "Amber"), label_visibility="collapsed")
        with col5:
            st.metric("🔴", sum(1 for d in st.session_state.decisions.values() if d == "Red"), label_visibility="collapsed")
        
        st.write("")
        
        if current < total:
            col_main, col_side = st.columns([2.5, 1])
            
            with col_main:
                smarts_pattern = df.iloc[current]['SMARTS']
                st.code(smarts_pattern, language='text')
                
                if 'Description' in df.columns and pd.notna(df.iloc[current]['Description']):
                    st.caption(df.iloc[current]['Description'])
                
                try:
                    pattern = Chem.MolFromSmarts(smarts_pattern)
                    if pattern is None:
                        st.error("⚠️ Invalid SMARTS")
                    else:
                        image_displayed = False
                        if use_smartsplus and st.session_state.api_key:
                            smartsplus_svg = get_smartsplus_image(smarts_pattern, st.session_state.api_key)
                            if smartsplus_svg:
                                # Remove width/height attributes so viewBox controls scaling
                                import re
                                svg_fixed = re.sub(r'\s*width="[^"]*"', '', smartsplus_svg)
                                svg_fixed = re.sub(r'\s*height="[^"]*"', '', svg_fixed)
                                
                                # Wrap in 800px container
                                scaled_svg = f"""
                                <div style="width: 800px; max-width: 100%;">
                                    {svg_fixed}
                                </div>
                                """
                                st.markdown(scaled_svg, unsafe_allow_html=True)
                                st.caption("🎨 SMARTS.plus")
                                image_displayed = True
                        
                        if not image_displayed:
                            img = Draw.MolToImage(pattern, size=(350, 280))
                            st.image(img, width=350)
                except Exception as e:
                    st.error(f"Error: {str(e)}")
            
            with col_side:
                if current in st.session_state.decisions:
                    color_map = {"OK": "🟢", "Yellow": "🟡", "Amber": "🟠", "Red": "🔴"}
                    st.info(f"{color_map.get(st.session_state.decisions[current], '')} {st.session_state.decisions[current]}")
                
                col_b1, col_b2 = st.columns(2)
                with col_b1:
                    if st.button("🟢 OK", use_container_width=True, key="ok"):
                        st.session_state.decisions[current] = "OK"
                        if current < total - 1:
                            st.session_state.current_idx += 1
                        st.rerun()
                    if st.button("🟠 Amber", use_container_width=True, key="amber"):
                        st.session_state.decisions[current] = "Amber"
                        if current < total - 1:
                            st.session_state.current_idx += 1
                        st.rerun()
                with col_b2:
                    if st.button("🟡 Yellow", use_container_width=True, key="yellow"):
                        st.session_state.decisions[current] = "Yellow"
                        if current < total - 1:
                            st.session_state.current_idx += 1
                        st.rerun()
                    if st.button("🔴 Red", use_container_width=True, key="red"):
                        st.session_state.decisions[current] = "Red"
                        if current < total - 1:
                            st.session_state.current_idx += 1
                        st.rerun()
                
                st.write("")
                
                col_n1, col_n2 = st.columns(2)
                with col_n1:
                    if st.button("⬅️", disabled=(current == 0), use_container_width=True):
                        st.session_state.current_idx -= 1
                        st.rerun()
                with col_n2:
                    if st.button("➡️", disabled=(current >= total - 1), use_container_width=True):
                        st.session_state.current_idx += 1
                        st.rerun()
                
                jump_to = st.number_input("Go to #", 1, total, current + 1, key='jump', label_visibility="collapsed")
                if st.button("Jump", use_container_width=True):
                    st.session_state.current_idx = jump_to - 1
                    st.rerun()
        
        else:
            st.success("🎉 Review complete!")
        
        if len(st.session_state.decisions) > 0:
            with st.expander("📥 Export Results"):
                results_df = df.copy()
                results_df['Decision'] = results_df.index.map(
                    lambda x: st.session_state.decisions.get(x, "NOT_REVIEWED")
                )
                
                col_e1, col_e2, col_e3 = st.columns(3)
                with col_e1:
                    st.download_button("📥 All", results_df.to_csv(index=False), "all_results.csv", "text/csv", use_container_width=True)
                with col_e2:
                    concerning = results_df[results_df['Decision'].isin(['Amber', 'Red'])]
                    if len(concerning) > 0:
                        st.download_button("⚠️ Concerning", concerning.to_csv(index=False), "concerning.csv", "text/csv", use_container_width=True)
                with col_e3:
                    ok_df = results_df[results_df['Decision'] == 'OK']
                    if len(ok_df) > 0:
                        st.download_button("✅ OK", ok_df.to_csv(index=False), "ok.csv", "text/csv", use_container_width=True)
    
    else:
        st.info("Upload a CSV with SMARTS patterns to begin")

# ============================================================================
# MODE 2: VALIDATOR
# ============================================================================

elif mode == "Validator":
    # Sub-mode selector
    validator_mode = st.radio(
        "Validation Type:",
        ["Batch Validation", "Quick Test"],
        horizontal=True,
        key="validator_submode"
    )
    
    st.write("")
    
    # ========================================================================
    # BATCH VALIDATION (Original Workflow)
    # ========================================================================
    if validator_mode == "Batch Validation":
        col_u1, col_u2 = st.columns(2)
        with col_u1:
            smarts_file = st.file_uploader("📁 SMARTS (CSV/SDF)", type=['csv', 'sdf'], key='val_smarts', label_visibility="collapsed")
        with col_u2:
            molecules_file = st.file_uploader("📁 Molecules (CSV/SDF)", type=['csv', 'sdf'], key='val_mols', label_visibility="collapsed")
        
        if smarts_file and st.session_state.smarts_data is None:
            try:
                if smarts_file.name.endswith('.csv'):
                    df = pd.read_csv(smarts_file)
                    if 'SMARTS' not in df.columns:
                        df.columns = ['SMARTS'] + list(df.columns[1:])
                else:
                    df = PandasTools.LoadSDF(smarts_file)
                st.session_state.smarts_data = df
                st.session_state.current_idx = 0
                st.session_state.decisions = {}
            except Exception as e:
                st.error(f"Error: {str(e)}")
        
        if molecules_file and st.session_state.test_molecules is None:
            try:
                if molecules_file.name.endswith('.csv'):
                    mol_df = pd.read_csv(molecules_file)
                    smiles_col = 'SMILES' if 'SMILES' in mol_df.columns else mol_df.columns[0]
                    mol_df['Mol'] = mol_df[smiles_col].apply(lambda x: Chem.MolFromSmiles(x) if pd.notna(x) else None)
                    mol_df = mol_df[mol_df['Mol'].notna()]
                else:
                    mol_df = PandasTools.LoadSDF(molecules_file)
                st.session_state.test_molecules = mol_df
                st.success(f"✅ Loaded {len(mol_df)} molecules")
            except Exception as e:
                st.error(f"Error: {str(e)}")
        
        if st.session_state.smarts_data is not None and st.session_state.test_molecules is not None:
            df = st.session_state.smarts_data
            mol_df = st.session_state.test_molecules
            total = len(df)
            current = st.session_state.current_idx
            
            st.progress(current / total if total > 0 else 0, text=f"Pattern {current + 1}/{total}")
            
            col_m1, col_m2 = st.columns(2)
            with col_m1:
                st.metric("Flagged", sum(1 for d in st.session_state.decisions.values() if d == "FLAGGED"))
            with col_m2:
                st.metric("Checked", sum(1 for d in st.session_state.decisions.values() if d == "CHECKED"))
            
            st.write("")
            
            if current < total:
                smarts_pattern = df.iloc[current]['SMARTS']
                
                col_main, col_side = st.columns([2.5, 1])
                
                with col_main:
                    st.code(smarts_pattern, language='text')
                    if 'Description' in df.columns and pd.notna(df.iloc[current]['Description']):
                        st.caption(df.iloc[current]['Description'])
                    
                    try:
                        pattern = Chem.MolFromSmarts(smarts_pattern)
                        if pattern:
                            matches = []
                            for idx, row in mol_df.iterrows():
                                mol = row.get('Mol') or row.get('ROMol')
                                if mol and mol.HasSubstructMatch(pattern):
                                    matches.append((idx, mol))
                            
                            st.write(f"**{len(matches)}/{len(mol_df)} matches ({len(matches)/len(mol_df)*100:.1f}%)**")
                            
                            if len(matches) > 0:
                                cols = st.columns(3)
                                for i, (idx, mol) in enumerate(matches[:6]):
                                    with cols[i % 3]:
                                        img = Draw.MolToImage(mol, size=(150, 150), 
                                                             highlightAtoms=mol.GetSubstructMatch(pattern))
                                        st.image(img, width=150)
                                        if 'Name' in mol_df.columns:
                                            st.caption(mol_df.iloc[idx]['Name'], unsafe_allow_html=True)
                            else:
                                st.info("No matches")
                    except Exception as e:
                        st.error(f"Error: {str(e)}")
                
                with col_side:
                    if current in st.session_state.decisions:
                        status = st.session_state.decisions[current]
                        if status == "FLAGGED":
                            st.warning("🚩 Flagged")
                        else:
                            st.info("✓ Checked")
                    
                    if st.button("🚩 FLAG", use_container_width=True, type="primary"):
                        st.session_state.decisions[current] = "FLAGGED"
                        st.rerun()
                    
                    if st.button("✓ OK", use_container_width=True):
                        st.session_state.decisions[current] = "CHECKED"
                        if current < total - 1:
                            st.session_state.current_idx += 1
                        st.rerun()
                    
                    st.write("")
                    col_n1, col_n2 = st.columns(2)
                    with col_n1:
                        if st.button("⬅️", disabled=(current == 0), use_container_width=True, key='prev_v'):
                            st.session_state.current_idx -= 1
                            st.rerun()
                    with col_n2:
                        if st.button("➡️", disabled=(current >= total - 1), use_container_width=True, key='next_v'):
                            st.session_state.current_idx += 1
                            st.rerun()
            
            else:
                st.success("🎉 Review complete!")
            
            if len(st.session_state.decisions) > 0:
                with st.expander("📥 Export"):
                    results_df = df.copy()
                    results_df['Status'] = results_df.index.map(lambda x: st.session_state.decisions.get(x, "NOT_REVIEWED"))
                    flagged = results_df[results_df['Status'] == 'FLAGGED']
                    
                    col_e1, col_e2 = st.columns(2)
                    with col_e1:
                        st.download_button("📥 All", results_df.to_csv(index=False), "validation.csv", "text/csv", use_container_width=True)
                    with col_e2:
                        if len(flagged) > 0:
                            st.download_button("🚩 Flagged", flagged.to_csv(index=False), "flagged.csv", "text/csv", use_container_width=True)
    
        else:
            st.info("Upload both SMARTS and molecules to begin")
    
    # ========================================================================
    # QUICK TEST (New Sub-Mode)
    # ========================================================================
    elif validator_mode == "Quick Test":
        # SMARTS input
        quick_smarts = st.text_area(
            "SMARTS Pattern:",
            height=80,
            placeholder="[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1",
            help="Paste your SMARTS pattern to test"
        )
        
        # Molecule source selector
        mol_source = st.radio(
            "Test Against:",
            ["Upload Molecule File", "Single SMILES"],
            horizontal=True,
            key="quick_mol_source"
        )
        
        st.write("")
        
        # Initialize session state for quick test
        if 'quick_test_results' not in st.session_state:
            st.session_state.quick_test_results = None
        if 'quick_page_idx' not in st.session_state:
            st.session_state.quick_page_idx = 0
        
        # Molecule input based on source
        if mol_source == "Upload Molecule File":
            quick_mol_file = st.file_uploader(
                "📁 Molecules (CSV/SDF)",
                type=['csv', 'sdf'],
                key='quick_mols',
                label_visibility="collapsed"
            )
            
            # Load molecules
            quick_mol_df = None
            if quick_mol_file:
                try:
                    if quick_mol_file.name.endswith('.csv'):
                        quick_mol_df = pd.read_csv(quick_mol_file)
                        smiles_col = 'SMILES' if 'SMILES' in quick_mol_df.columns else quick_mol_df.columns[0]
                        quick_mol_df['Mol'] = quick_mol_df[smiles_col].apply(
                            lambda x: Chem.MolFromSmiles(x) if pd.notna(x) else None
                        )
                        quick_mol_df = quick_mol_df[quick_mol_df['Mol'].notna()]
                    else:
                        quick_mol_df = PandasTools.LoadSDF(quick_mol_file, molColName='Mol')
                    
                    st.success(f"✅ Loaded {len(quick_mol_df)} molecules")
                except Exception as e:
                    st.error(f"Error loading molecules: {str(e)}")
            
            # Run test button
            if st.button("🚀 Run Test", type="primary", use_container_width=True, disabled=not (quick_smarts and quick_mol_df is not None)):
                if quick_smarts and quick_mol_df is not None:
                    with st.spinner("Testing SMARTS pattern..."):
                        try:
                            pattern = Chem.MolFromSmarts(quick_smarts)
                            if pattern is None:
                                st.error("❌ Invalid SMARTS pattern")
                            else:
                                # Find matches
                                matches = []
                                for idx, row in quick_mol_df.iterrows():
                                    mol = row.get('Mol') or row.get('ROMol')
                                    if mol and mol.HasSubstructMatch(pattern):
                                        matches.append((idx, mol, row))
                                
                                st.session_state.quick_test_results = {
                                    'pattern': pattern,
                                    'smarts': quick_smarts,
                                    'matches': matches,
                                    'total': len(quick_mol_df),
                                    'mol_df': quick_mol_df
                                }
                                st.session_state.quick_page_idx = 0
                                st.rerun()
                        except Exception as e:
                            st.error(f"Error: {str(e)}")
        
        else:  # Single SMILES mode
            quick_smiles = st.text_input(
                "SMILES:",
                placeholder="CCO",
                help="Paste a single SMILES string"
            )
            
            # Run test button
            if st.button("🚀 Test Match", type="primary", use_container_width=True, disabled=not (quick_smarts and quick_smiles)):
                if quick_smarts and quick_smiles:
                    try:
                        pattern = Chem.MolFromSmarts(quick_smarts)
                        mol = Chem.MolFromSmiles(quick_smiles)
                        
                        if pattern is None:
                            st.error("❌ Invalid SMARTS pattern")
                        elif mol is None:
                            st.error("❌ Invalid SMILES")
                        else:
                            is_match = mol.HasSubstructMatch(pattern)
                            
                            col_res1, col_res2 = st.columns([1, 2])
                            
                            with col_res1:
                                if is_match:
                                    st.success("✅ MATCH")
                                else:
                                    st.error("❌ NO MATCH")
                            
                            with col_res2:
                                if is_match:
                                    match_atoms = mol.GetSubstructMatch(pattern)
                                    img = Draw.MolToImage(mol, size=(300, 300), highlightAtoms=match_atoms)
                                else:
                                    img = Draw.MolToImage(mol, size=(300, 300))
                                st.image(img, width=300)
                    
                    except Exception as e:
                        st.error(f"Error: {str(e)}")
        
        # Display results for file mode
        if st.session_state.quick_test_results and mol_source == "Upload Molecule File":
            results = st.session_state.quick_test_results
            matches = results['matches']
            total = results['total']
            
            st.write("---")
            st.subheader("📊 Results")
            
            # Summary
            col_r1, col_r2, col_r3 = st.columns(3)
            with col_r1:
                st.metric("Total Molecules", total)
            with col_r2:
                st.metric("✅ Matches", len(matches))
            with col_r3:
                st.metric("Match Rate", f"{len(matches)/total*100:.1f}%")
            
            st.write("")
            
            if len(matches) > 0:
                # Pagination
                per_page = 6
                total_pages = (len(matches) + per_page - 1) // per_page
                current_page = st.session_state.quick_page_idx
                
                # Display matches
                start_idx = current_page * per_page
                end_idx = min(start_idx + per_page, len(matches))
                page_matches = matches[start_idx:end_idx]
                
                cols = st.columns(3)
                for i, (idx, mol, row) in enumerate(page_matches):
                    with cols[i % 3]:
                        match_atoms = mol.GetSubstructMatch(results['pattern'])
                        img = Draw.MolToImage(mol, size=(150, 150), highlightAtoms=match_atoms)
                        st.image(img, width=150)
                        
                        # Show molecule name if available
                        if 'Name' in row:
                            st.caption(row['Name'])
                        elif 'SMILES' in row:
                            smiles_short = row['SMILES'][:30] + "..." if len(row['SMILES']) > 30 else row['SMILES']
                            st.caption(smiles_short)
                
                # Pagination controls
                if total_pages > 1:
                    st.write("")
                    col_p1, col_p2, col_p3 = st.columns([1, 2, 1])
                    
                    with col_p1:
                        if st.button("⬅️ Previous", disabled=(current_page == 0), use_container_width=True):
                            st.session_state.quick_page_idx -= 1
                            st.rerun()
                    
                    with col_p2:
                        st.write(f"Page {current_page + 1} of {total_pages}")
                    
                    with col_p3:
                        if st.button("Next ➡️", disabled=(current_page >= total_pages - 1), use_container_width=True):
                            st.session_state.quick_page_idx += 1
                            st.rerun()
                
                # Export matches
                st.write("")
                if st.button("📥 Export Matches", use_container_width=True):
                    # Create export dataframe
                    match_indices = [idx for idx, _, _ in matches]
                    export_df = results['mol_df'].iloc[match_indices].copy()
                    
                    # Remove Mol column for CSV export
                    export_cols = [col for col in export_df.columns if col not in ['Mol', 'ROMol']]
                    
                    # Add SMILES if not present
                    if 'SMILES' not in export_cols:
                        export_df['SMILES'] = export_df['Mol'].apply(lambda m: Chem.MolToSmiles(m) if m else '')
                        export_cols = ['SMILES'] + export_cols
                    
                    st.download_button(
                        "📥 Download Matches CSV",
                        export_df[export_cols].to_csv(index=False),
                        "smarts_matches.csv",
                        "text/csv",
                        use_container_width=True
                    )
            else:
                st.info("No matches found for this SMARTS pattern")


# ============================================================================
# MODE 3: GEN AI FILTER (NEW!)
# ============================================================================

elif mode == "Gen AI Filter":
    st.subheader("🧬 Batch Filter Gen AI Output")
    
    col_u1, col_u2 = st.columns(2)
    with col_u1:
        filter_smarts_file = st.file_uploader(
            "📁 Approved SMARTS (CSV from Visualizer)",
            type=['csv'],
            key='filter_smarts',
            help="Upload CSV with 'Decision' column (from Visualizer export)"
        )
    with col_u2:
        gen_ai_file = st.file_uploader(
            "📁 Gen AI Molecules (SDF or CSV)",
            type=['sdf', 'csv'],
            key='gen_ai_mols',
            help="SDF with molecules or CSV with SMILES column"
        )
    
    # Filter severity selector
    filter_severity = st.multiselect(
        "Apply patterns with these decisions:",
        ["Red", "Amber", "Yellow", "OK"],
        default=["Red", "Amber"],
        help="Which traffic light categories to use as filters"
    )
    
    # Load filter patterns
    filter_patterns = None
    if filter_smarts_file:
        try:
            filter_df = pd.read_csv(filter_smarts_file)
            
            # Check if Decision column exists
            if 'Decision' not in filter_df.columns:
                st.warning("No 'Decision' column found. Using all patterns.")
                filter_patterns = filter_df
            else:
                # Filter by selected severities
                filter_patterns = filter_df[filter_df['Decision'].isin(filter_severity)]
                st.info(f"✅ Loaded {len(filter_patterns)} filter patterns ({', '.join(filter_severity)})")
                
        except Exception as e:
            st.error(f"Error loading SMARTS: {str(e)}")
    
    # Load molecules
    gen_ai_mols = None
    if gen_ai_file:
        try:
            if gen_ai_file.name.endswith('.sdf'):
                # Load SDF
                gen_ai_mols = PandasTools.LoadSDF(gen_ai_file, molColName='Mol')
                st.success(f"✅ Loaded {len(gen_ai_mols)} molecules from SDF")
            else:
                # Load CSV with SMILES
                gen_ai_mols = pd.read_csv(gen_ai_file)
                smiles_col = 'SMILES' if 'SMILES' in gen_ai_mols.columns else gen_ai_mols.columns[0]
                gen_ai_mols['Mol'] = gen_ai_mols[smiles_col].apply(
                    lambda x: Chem.MolFromSmiles(x) if pd.notna(x) else None
                )
                gen_ai_mols = gen_ai_mols[gen_ai_mols['Mol'].notna()]
                st.success(f"✅ Loaded {len(gen_ai_mols)} molecules from CSV")
                
        except Exception as e:
            st.error(f"Error loading molecules: {str(e)}")
    
    # Run filtering
    if filter_patterns is not None and gen_ai_mols is not None and len(filter_patterns) > 0:
        
        if st.button("🚀 Run Batch Filter", type="primary", use_container_width=True):
            
            with st.spinner("Filtering molecules..."):
                
                # Compile SMARTS patterns
                compiled_patterns = []
                for idx, row in filter_patterns.iterrows():
                    smarts = row['SMARTS']
                    pattern = Chem.MolFromSmarts(smarts)
                    if pattern:
                        compiled_patterns.append({
                            'pattern': pattern,
                            'smarts': smarts,
                            'description': row.get('Description', ''),
                            'decision': row.get('Decision', 'Unknown')
                        })
                
                # Filter molecules
                passed = []
                failed = []
                rejection_reasons = []
                
                for idx, row in gen_ai_mols.iterrows():
                    mol = row.get('Mol') or row.get('ROMol')
                    if mol is None:
                        continue
                    
                    caught_by = []
                    for p in compiled_patterns:
                        if mol.HasSubstructMatch(p['pattern']):
                            caught_by.append(p)
                    
                    if len(caught_by) == 0:
                        passed.append(idx)
                    else:
                        failed.append(idx)
                        # Get SMILES for the rejected molecule
                        mol_smiles = Chem.MolToSmiles(mol) if mol else ''
                        rejection_reasons.append({
                            'molecule_idx': idx,
                            'SMILES': mol_smiles,
                            'num_violations': len(caught_by),
                            'patterns': [p['smarts'] for p in caught_by],
                            'descriptions': [p['description'] for p in caught_by]
                        })
                
                # Store results
                st.session_state.filter_results = {
                    'passed': passed,
                    'failed': failed,
                    'rejection_reasons': rejection_reasons,
                    'total': len(gen_ai_mols),
                    'patterns_used': len(compiled_patterns)
                }
    
    # Display results
    if st.session_state.filter_results:
        results = st.session_state.filter_results
        
        st.write("---")
        st.subheader("📊 Filtering Results")
        
        # Summary metrics
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Total Molecules", results['total'])
        with col2:
            st.metric("✅ Passed", len(results['passed']), 
                     delta=f"{len(results['passed'])/results['total']*100:.1f}%")
        with col3:
            st.metric("❌ Failed", len(results['failed']),
                     delta=f"{len(results['failed'])/results['total']*100:.1f}%")
        with col4:
            st.metric("Filters Used", results['patterns_used'])
        
        # Export buttons
        st.write("")
        col_e1, col_e2 = st.columns(2)
        
        with col_e1:
            # Export passed molecules
            if len(results['passed']) > 0:
                passed_mols = gen_ai_mols.iloc[results['passed']].copy()
                
                # Create export dataframe
                if 'SMILES' in passed_mols.columns:
                    export_df = passed_mols[['SMILES']].copy()
                else:
                    # Generate SMILES from Mol objects
                    export_df = pd.DataFrame({
                        'SMILES': passed_mols['Mol'].apply(lambda m: Chem.MolToSmiles(m) if m else '')
                    })
                
                # Add other columns if they exist
                for col in gen_ai_mols.columns:
                    if col not in ['Mol', 'ROMol', 'SMILES'] and col in passed_mols.columns:
                        export_df[col] = passed_mols[col].values
                
                st.download_button(
                    "📥 Download Passed Molecules",
                    export_df.to_csv(index=False),
                    "passed_molecules.csv",
                    "text/csv",
                    use_container_width=True,
                    type="primary"
                )
        
        with col_e2:
            # Export rejection report
            if len(results['failed']) > 0:
                rejection_df = pd.DataFrame(results['rejection_reasons'])
                
                st.download_button(
                    "📥 Download Rejection Report",
                    rejection_df.to_csv(index=False),
                    "rejection_report.csv",
                    "text/csv",
                    use_container_width=True
                )
        
        # Show sample rejections
        if len(results['failed']) > 0:
            with st.expander(f"🔍 View Sample Rejections (first 5 of {len(results['failed'])})"):
                for i, reason in enumerate(results['rejection_reasons'][:5]):
                    st.write(f"**Molecule {reason['molecule_idx']}** - Caught by {reason['num_violations']} pattern(s):")
                    for j, (pattern, desc) in enumerate(zip(reason['patterns'], reason['descriptions'])):
                        st.text(f"  • {pattern}")
                        if desc:
                            st.caption(f"    {desc}")
                    if i < 4:
                        st.divider()
    
    elif filter_patterns is not None and gen_ai_mols is not None and len(filter_patterns) == 0:
        st.warning("No patterns selected. Choose at least one decision category to filter.")
    
    else:
        st.info("Upload both approved SMARTS patterns and Gen AI molecules, then click 'Run Batch Filter'")
