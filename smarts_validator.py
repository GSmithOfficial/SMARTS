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

# NEW: Quick Guide expander
with st.expander("❓ Quick Guide & Resources", expanded=False):
    col_g1, col_g2 = st.columns(2)
    
    with col_g1:
        st.markdown("""
        ### 📖 SMARTS Resources
        **Learn SMARTS Pattern Language:**
        - [Daylight SMARTS Tutorial](https://www.daylight.com/dayhtml_tutorials/languages/smarts/) - Official reference
        - [SMARTS Examples](https://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html) - Common patterns
        - [SMARTS Theory Manual](https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html) - Deep dive
        
        ### 🔬 Mode Overview
        **Visualizer** - Review and curate SMARTS pattern libraries  
        **Validator** - Test pattern specificity against molecules  
        **Gen AI Filter** - Batch filter molecules with curated patterns
        """)
    
    with col_g2:
        st.markdown("""
        ### 📁 Expected File Formats
        **SMARTS CSV**: Must contain `SMARTS` column  
        - Optional: `Description`, `Decision` columns
        
        **Molecule Files**: 
        - `.sdf` format (direct structure import)
        - `.csv` with `SMILES` column
        
        **Filter Patterns**: CSV with `SMARTS` + `Decision` columns
        
        ### 🚦 Traffic Light System
        🔴 **Red** - Critical liabilities (always block)  
        🟠 **Amber** - Cautionary flags (context-dependent)  
        🟡 **Yellow** - Minor concerns (review recommended)  
        🟢 **OK** - Acceptable patterns
        """)

# Mode selector
mode = st.radio("Select Mode:", ["Visualizer", "Validator", "Gen AI Filter"], horizontal=True, key="mode_selector")

# NEW: Add contextual description under mode selector
mode_descriptions = {
    "Visualizer": "Review and classify SMARTS patterns from your library",
    "Validator": "Test SMARTS patterns against molecule sets to check specificity",
    "Gen AI Filter": "Filter generated molecules using curated SMARTS libraries"
}
st.caption(mode_descriptions[mode])

st.session_state.mode = mode

# Settings expander - UPDATED with disclaimer
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
            disabled=not bool(st.session_state.api_key),
            help="Enable high-quality SMARTS visualizations"
        )
        st.session_state.debug_mode = st.checkbox(
            "Debug mode",
            value=False,
            help="Show SMARTS.plus API errors"
        )
    
    # NEW: API Credits & Terms section
    st.divider()
    st.markdown("#### API Credits & Terms")
    st.info("""
**Visualization powered by [SMARTS.plus](https://smarts.plus)**  
• Free API access requires registration at [smarts.plus/sign_up](https://smarts.plus/sign_up)  
• By using SMARTS.plus visualization, you agree to their [usage policies](https://api.smarts.plus/usage_policies/)  
• API keys are personal and non-transferable  
• External applications require API key for POST requests

*This tool is independent software - not affiliated with or endorsed by SMARTS.plus.*
    """)

st.divider()

# ============================================================================
# MODE 1: VISUALIZER
# ============================================================================

if mode == "Visualizer":
    # Create tabs for different input methods
    tab1, tab2 = st.tabs(["📁 CSV Upload", "✍️ Single SMARTS"])
    
    # ========================================================================
    # TAB 1: CSV UPLOAD (Original functionality)
    # ========================================================================
    with tab1:
        uploaded_file = st.file_uploader("Upload SMARTS CSV", type=['csv'], label_visibility="collapsed")
        
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
                st.metric("Total", total)
            with col2:
                st.metric("Reviewed", len(st.session_state.decisions))
            with col3:
                st.metric("Remaining", total - len(st.session_state.decisions))
            with col4:
                red_count = sum(1 for d in st.session_state.decisions.values() if d == "Red")
                st.metric("🔴 Red", red_count)
            with col5:
                amber_count = sum(1 for d in st.session_state.decisions.values() if d == "Amber")
                st.metric("🟠 Amber", amber_count)
            
            st.write("")
            
            # Navigation
            col_nav1, col_nav2, col_nav3, col_nav4 = st.columns([1, 1, 1, 1])
            with col_nav1:
                if st.button("⏮️ First", use_container_width=True):
                    st.session_state.current_idx = 0
                    st.rerun()
            with col_nav2:
                if st.button("◀️ Previous", use_container_width=True):
                    if current > 0:
                        st.session_state.current_idx = current - 1
                        st.rerun()
            with col_nav3:
                if st.button("▶️ Next", use_container_width=True):
                    if current < total - 1:
                        st.session_state.current_idx = current + 1
                        st.rerun()
            with col_nav4:
                if st.button("⏭️ Last", use_container_width=True):
                    st.session_state.current_idx = total - 1
                    st.rerun()
            
            st.write("")
            
            # Current pattern display
            if 0 <= current < total:
                row = df.iloc[current]
                smarts = row['SMARTS']
                
                col_info1, col_info2 = st.columns([1, 2])
                with col_info1:
                    st.text_input("SMARTS Pattern", value=smarts, disabled=True, key="current_smarts")
                with col_info2:
                    if 'Description' in row and pd.notna(row['Description']):
                        st.text_input("Description", value=row['Description'], disabled=True)
                
                # Visualization
                st.write("")
                
                # Try SMARTS.plus first if enabled
                if use_smartsplus and st.session_state.api_key:
                    svg_data = get_smartsplus_image(smarts, st.session_state.api_key)
                    if svg_data:
                        st.markdown(svg_data, unsafe_allow_html=True)
                    else:
                        # Fallback to RDKit
                        mol = Chem.MolFromSmarts(smarts)
                        if mol:
                            img = Draw.MolToImage(mol, size=(400, 200))
                            st.image(img, use_container_width=False)
                        else:
                            st.error("Invalid SMARTS pattern")
                else:
                    # Use RDKit visualization
                    mol = Chem.MolFromSmarts(smarts)
                    if mol:
                        img = Draw.MolToImage(mol, size=(400, 200))
                        st.image(img, use_container_width=False)
                    else:
                        st.error("Invalid SMARTS pattern")
                
                st.write("")
                
                # Decision buttons
                st.write("**Classify this pattern:**")
                col_btn1, col_btn2, col_btn3, col_btn4 = st.columns(4)
                
                current_decision = st.session_state.decisions.get(current, None)
                
                with col_btn1:
                    red_type = "primary" if current_decision == "Red" else "secondary"
                    if st.button("🔴 Red", type=red_type, use_container_width=True):
                        st.session_state.decisions[current] = "Red"
                        if current < total - 1:
                            st.session_state.current_idx = current + 1
                        st.rerun()
                
                with col_btn2:
                    amber_type = "primary" if current_decision == "Amber" else "secondary"
                    if st.button("🟠 Amber", type=amber_type, use_container_width=True):
                        st.session_state.decisions[current] = "Amber"
                        if current < total - 1:
                            st.session_state.current_idx = current + 1
                        st.rerun()
                
                with col_btn3:
                    yellow_type = "primary" if current_decision == "Yellow" else "secondary"
                    if st.button("🟡 Yellow", type=yellow_type, use_container_width=True):
                        st.session_state.decisions[current] = "Yellow"
                        if current < total - 1:
                            st.session_state.current_idx = current + 1
                        st.rerun()
                
                with col_btn4:
                    ok_type = "primary" if current_decision == "OK" else "secondary"
                    if st.button("🟢 OK", type=ok_type, use_container_width=True):
                        st.session_state.decisions[current] = "OK"
                        if current < total - 1:
                            st.session_state.current_idx = current + 1
                        st.rerun()
                
                if current_decision:
                    st.success(f"Current decision: {current_decision}")
                
                # Export results
                if len(st.session_state.decisions) > 0:
                    st.write("")
                    st.write("---")
                    
                    export_df = df.copy()
                    export_df['Decision'] = export_df.index.map(
                        lambda i: st.session_state.decisions.get(i, "")
                    )
                    
                    csv = export_df.to_csv(index=False)
                    st.download_button(
                        "📥 Download Classified Patterns",
                        csv,
                        "classified_smarts.csv",
                        "text/csv",
                        use_container_width=True,
                        type="primary"
                    )
    
    # ========================================================================
    # TAB 2: SINGLE SMARTS (Quick check)
    # ========================================================================
    with tab2:
        st.write("**Quick SMARTS pattern check**")
        
        single_smarts = st.text_input("Enter SMARTS pattern:", key="single_smarts_input")
        
        if single_smarts:
            col_v1, col_v2 = st.columns([1, 2])
            
            with col_v1:
                mol = Chem.MolFromSmarts(single_smarts)
                if mol:
                    st.success("✅ Valid SMARTS")
                    st.write(f"**Atoms:** {mol.GetNumAtoms()}")
                    st.write(f"**Bonds:** {mol.GetNumBonds()}")
                else:
                    st.error("❌ Invalid SMARTS")
            
            with col_v2:
                if mol:
                    # Try SMARTS.plus first if enabled
                    if use_smartsplus and st.session_state.api_key:
                        svg_data = get_smartsplus_image(single_smarts, st.session_state.api_key, use_cache=False)
                        if svg_data:
                            st.markdown(svg_data, unsafe_allow_html=True)
                        else:
                            img = Draw.MolToImage(mol, size=(400, 200))
                            st.image(img, use_container_width=False)
                    else:
                        img = Draw.MolToImage(mol, size=(400, 200))
                        st.image(img, use_container_width=False)

# ============================================================================
# MODE 2: VALIDATOR
# ============================================================================

elif mode == "Validator":
    st.write("**Test SMARTS patterns against molecules**")
    
    col_val1, col_val2 = st.columns(2)
    
    with col_val1:
        test_smarts = st.text_area("SMARTS Pattern:", height=100, key="validator_smarts")
    
    with col_val2:
        test_file = st.file_uploader("Molecules (CSV/SDF)", type=['csv', 'sdf'], key="validator_molecules")
    
    if test_smarts and test_file:
        # Parse SMARTS
        pattern = Chem.MolFromSmarts(test_smarts)
        
        if not pattern:
            st.error("❌ Invalid SMARTS pattern")
            st.stop()
        
        st.success("✅ Valid SMARTS pattern")
        
        # Load molecules
        try:
            if test_file.name.endswith('.sdf'):
                test_mols = PandasTools.LoadSDF(test_file, molColName='Mol')
            else:
                test_mols = pd.read_csv(test_file)
                smiles_col = 'SMILES' if 'SMILES' in test_mols.columns else test_mols.columns[0]
                test_mols['Mol'] = test_mols[smiles_col].apply(
                    lambda x: Chem.MolFromSmiles(x) if pd.notna(x) else None
                )
                test_mols = test_mols[test_mols['Mol'].notna()]
            
            st.success(f"✅ Loaded {len(test_mols)} molecules")
        except Exception as e:
            st.error(f"Error loading molecules: {str(e)}")
            st.stop()
        
        # Run matching
        if st.button("🔍 Test Pattern", type="primary", use_container_width=True):
            matches = []
            non_matches = []
            
            with st.spinner("Testing pattern..."):
                for idx, row in test_mols.iterrows():
                    mol = row.get('Mol') or row.get('ROMol')
                    if mol and mol.HasSubstructMatch(pattern):
                        matches.append(idx)
                    else:
                        non_matches.append(idx)
            
            # Results
            st.write("---")
            col_r1, col_r2, col_r3 = st.columns(3)
            with col_r1:
                st.metric("Total Tested", len(test_mols))
            with col_r2:
                match_pct = f"{len(matches)/len(test_mols)*100:.1f}%"
                st.metric("✅ Matches", len(matches), delta=match_pct)
            with col_r3:
                st.metric("❌ No Match", len(non_matches))
            
            # Show sample matches
            if len(matches) > 0:
                with st.expander(f"🔍 View Sample Matches (first 5 of {len(matches)})"):
                    sample_matches = test_mols.iloc[matches[:5]]
                    for idx, row in sample_matches.iterrows():
                        mol = row.get('Mol') or row.get('ROMol')
                        if mol:
                            img = Draw.MolToImage(mol, size=(300, 150), highlightAtoms=mol.GetSubstructMatch(pattern))
                            st.image(img, caption=f"Molecule {idx}")
                            st.divider()

# ============================================================================
# MODE 3: GEN AI FILTER
# ============================================================================

elif mode == "Gen AI Filter":
    st.write("**Filter AI-generated molecules with curated SMARTS patterns**")
    
    col_filter1, col_filter2 = st.columns(2)
    
    with col_filter1:
        filter_smarts_file = st.file_uploader(
            "Approved SMARTS Patterns (CSV with Decision column)",
            type=['csv'],
            key="filter_smarts"
        )
    
    with col_filter2:
        gen_ai_file = st.file_uploader(
            "Gen AI Molecules (CSV/SDF)",
            type=['csv', 'sdf'],
            key="gen_ai_mols"
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
                # Load CSV with SMILES - handle Excel's sep= line
                # Check if first line is "sep=,"
                gen_ai_file.seek(0)
                first_line = gen_ai_file.readline().decode('utf-8').strip()
                gen_ai_file.seek(0)
                
                if first_line.startswith('sep='):
                    # Skip the sep= line
                    gen_ai_mols = pd.read_csv(gen_ai_file, skiprows=1)
                else:
                    gen_ai_mols = pd.read_csv(gen_ai_file)
                
                # Find SMILES column
                smiles_col = 'SMILES' if 'SMILES' in gen_ai_mols.columns else gen_ai_mols.columns[0]
                
                # Convert SMILES to molecules
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
                            'descriptions': [p['description'] for p in caught_by],
                            'decisions': [p['decision'] for p in caught_by]
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
            if results['total'] > 0:
                pass_pct = f"{len(results['passed'])/results['total']*100:.1f}%"
            else:
                pass_pct = "N/A"
            st.metric("✅ Passed", len(results['passed']), delta=pass_pct)
        with col3:
            if results['total'] > 0:
                fail_pct = f"{len(results['failed'])/results['total']*100:.1f}%"
            else:
                fail_pct = "N/A"
            st.metric("❌ Failed", len(results['failed']), delta=fail_pct)
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
                    for j, (pattern, desc, decision) in enumerate(zip(reason['patterns'], reason['descriptions'], reason['decisions'])):
                        st.text(f"  • [{decision}] {pattern}")
                        if desc:
                            st.caption(f"    {desc}")
                    if i < 4:
                        st.divider()
    
    elif filter_patterns is not None and gen_ai_mols is not None and len(filter_patterns) == 0:
        st.warning("No patterns selected. Choose at least one decision category to filter.")
    
    else:
        st.info("Upload both approved SMARTS patterns and Gen AI molecules, then click 'Run Batch Filter'")
