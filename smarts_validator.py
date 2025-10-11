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
    /* Reduce padding and margins */
    .block-container {
        padding-top: 1rem;
        padding-bottom: 0rem;
        padding-left: 2rem;
        padding-right: 2rem;
    }
    
    /* Compact metrics */
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
    
    /* Compact buttons */
    .stButton button {
        padding: 0.25rem 0.5rem;
        font-size: 0.85rem;
    }
    
    /* Smaller headers */
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
    
    /* Compact file uploader */
    [data-testid="stFileUploader"] {
        padding: 0.5rem;
    }
    
    /* Reduce spacing between elements */
    .element-container {
        margin-bottom: 0.5rem;
    }
    
    /* Compact sidebar */
    .css-1d391kg {
        padding: 1rem 0.5rem;
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

# Helper function for SMARTS.plus visualization
def get_smartsplus_image(smarts_pattern, api_key, use_cache=True):
    """Get visualization from SMARTS.plus API"""
    cache_key = f"{smarts_pattern}_{api_key[:8]}"
    if use_cache and cache_key in st.session_state.viz_cache:
        return st.session_state.viz_cache[cache_key]
    
    try:
        payload = {
            "query": {
                "smarts": smarts_pattern,
                "parameters": {
                    "file_format": "png",
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
                        image_data = base64.b64decode(result['result']['image'])
                        st.session_state.viz_cache[cache_key] = image_data
                        return image_data
        elif response.status_code == 200:
            result = response.json()
            if 'result' in result and 'image' in result['result']:
                image_data = base64.b64decode(result['result']['image'])
                st.session_state.viz_cache[cache_key] = image_data
                return image_data
        return None
    except:
        return None

# Compact header
col_title, col_mode = st.columns([3, 2])
with col_title:
    st.title("🔬 SMARTS Toolkit")
with col_mode:
    mode = st.radio("Mode:", ["Visualizer", "Validator"], horizontal=True, label_visibility="collapsed")
    st.session_state.mode = mode

# Settings expander (collapsed by default)
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

st.divider()

# ============================================================================
# MODE 1: VISUALIZER
# ============================================================================

if mode == "Visualizer":
    # Compact file upload
    uploaded_file = st.file_uploader("📁 Upload SMARTS CSV", type=['csv'], label_visibility="collapsed")
    
    if uploaded_file:
        if st.session_state.smarts_data is None:
            try:
                df = pd.read_csv(uploaded_file)
                if 'SMARTS' not in df.columns:
                    df.columns = ['SMARTS'] + list(df.columns[1:])
                st.session_state.smarts_data = df
                st.session_state.current_idx = 0
                st.session_state.decisions = {}
            except Exception as e:
                st.error(f"Parse error: {str(e)}")
                st.stop()
        
        df = st.session_state.smarts_data
        total = len(df)
        current = st.session_state.current_idx
        progress = len(st.session_state.decisions) / total if total > 0 else 0
        
        # Compact progress bar
        st.progress(progress, text=f"Pattern {current + 1}/{total}")
        
        # Compact metrics
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
        
        st.write("")  # Small spacer
        
        if current < total:
            # Main layout - tighter columns
            col_main, col_side = st.columns([2.5, 1])
            
            with col_main:
                smarts_pattern = df.iloc[current]['SMARTS']
                
                # Pattern display
                st.code(smarts_pattern, language='text')
                
                if 'Description' in df.columns and pd.notna(df.iloc[current]['Description']):
                    st.caption(df.iloc[current]['Description'])
                
                # Visualization
                try:
                    pattern = Chem.MolFromSmarts(smarts_pattern)
                    if pattern is None:
                        st.error("⚠️ Invalid SMARTS")
                    else:
                        # Try SMARTS.plus
                        image_displayed = False
                        if use_smartsplus and st.session_state.api_key:
                            smartsplus_image = get_smartsplus_image(smarts_pattern, st.session_state.api_key)
                            if smartsplus_image:
                                st.image(smartsplus_image, width=400)
                                st.caption("🎨 SMARTS.plus")
                                image_displayed = True
                        
                        # Fallback to RDKit
                        if not image_displayed:
                            img = Draw.MolToImage(pattern, size=(350, 280))
                            st.image(img, width=350)
                except Exception as e:
                    st.error(f"Error: {str(e)}")
            
            with col_side:
                # Current status
                if current in st.session_state.decisions:
                    color_map = {"OK": "🟢", "Yellow": "🟡", "Amber": "🟠", "Red": "🔴"}
                    st.info(f"{color_map.get(st.session_state.decisions[current], '')} {st.session_state.decisions[current]}")
                
                # Decision buttons - 2x2 grid
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
                
                # Navigation
                col_n1, col_n2 = st.columns(2)
                with col_n1:
                    if st.button("⬅️", disabled=(current == 0), use_container_width=True):
                        st.session_state.current_idx -= 1
                        st.rerun()
                with col_n2:
                    if st.button("➡️", disabled=(current >= total - 1), use_container_width=True):
                        st.session_state.current_idx += 1
                        st.rerun()
                
                # Jump to
                jump_to = st.number_input("Go to #", 1, total, current + 1, key='jump', label_visibility="collapsed")
                if st.button("Jump", use_container_width=True):
                    st.session_state.current_idx = jump_to - 1
                    st.rerun()
        
        else:
            st.success("🎉 Review complete!")
        
        # Export (collapsible)
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
    col_u1, col_u2 = st.columns(2)
    with col_u1:
        smarts_file = st.file_uploader("📁 SMARTS (CSV/SDF)", type=['csv', 'sdf'], key='val_smarts', label_visibility="collapsed")
    with col_u2:
        molecules_file = st.file_uploader("📁 Molecules (CSV/SDF)", type=['csv', 'sdf'], key='val_mols', label_visibility="collapsed")
    
    # Load files
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
