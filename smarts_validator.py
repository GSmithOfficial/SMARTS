import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem import Descriptors
import io

# Page config
st.set_page_config(page_title="SMARTS Validator", layout="wide")

# Initialize session state
if 'current_idx' not in st.session_state:
    st.session_state.current_idx = 0
if 'decisions' not in st.session_state:
    st.session_state.decisions = {}
if 'smarts_data' not in st.session_state:
    st.session_state.smarts_data = None

# File upload
st.title("🔬 SMARTS Pattern Validator")

uploaded_file = st.file_uploader("Upload your SMARTS list (CSV with 'SMARTS' and optional 'Description' columns)", 
                                  type=['csv', 'txt'])

if uploaded_file:
    # Load SMARTS data
    if st.session_state.smarts_data is None:
        try:
            df = pd.read_csv(uploaded_file)
            if 'SMARTS' not in df.columns:
                # Assume first column is SMARTS if no header
                df.columns = ['SMARTS'] + list(df.columns[1:])
            st.session_state.smarts_data = df
        except:
            st.error("Could not parse file. Ensure CSV has 'SMARTS' column")
            st.stop()
    
    df = st.session_state.smarts_data
    total = len(df)
    current = st.session_state.current_idx
    
    # Progress bar
    progress = len(st.session_state.decisions) / total
    st.progress(progress)
    st.write(f"**Progress:** {len(st.session_state.decisions)}/{total} reviewed "
             f"({int(progress*100)}%)")
    
    # Display current SMARTS
    if current < total:
        col1, col2 = st.columns([2, 1])
        
        with col1:
            st.subheader(f"Pattern {current + 1} of {total}")
            smarts_pattern = df.iloc[current]['SMARTS']
            st.code(smarts_pattern, language='text')
            
            if 'Description' in df.columns:
                st.write(f"**Description:** {df.iloc[current]['Description']}")
            
            # Try to parse SMARTS
            try:
                pattern = Chem.MolFromSmarts(smarts_pattern)
                if pattern is None:
                    st.error("⚠️ INVALID SMARTS - Cannot parse!")
                else:
                    st.success("✅ Valid SMARTS syntax")
                    
                    # Show example molecules (you'd replace with your library)
                    st.write("**Example matches from test set:**")
                    
                    # Generate or load test molecules
                    test_smiles = [
                        "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O",  # Ibuprofen
                        "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
                        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # Caffeine
                        "c1ccc2c(c1)ccc3c2cccc3",  # Anthracene
                        "C1=CC=C(C=C1)C(=O)O",  # Benzoic acid
                        "CCO",  # Ethanol
                        "c1ccccc1",  # Benzene
                        "CC(C)(C)c1ccc(O)cc1"  # BHT fragment
                    ]
                    
                    matches_found = 0
                    mol_images = []
                    
                    for smi in test_smiles:
                        mol = Chem.MolFromSmiles(smi)
                        if mol and mol.HasSubstructMatch(pattern):
                            matches_found += 1
                            mol_images.append(Draw.MolToImage(mol, highlightAtoms=mol.GetSubstructMatch(pattern)))
                            if matches_found >= 5:  # Show max 5 examples
                                break
                    
                    if matches_found > 0:
                        cols = st.columns(min(matches_found, 3))
                        for i, img in enumerate(mol_images[:3]):
                            with cols[i]:
                                st.image(img, width=200)
                        st.write(f"*Matched {matches_found} molecules in test set*")
                    else:
                        st.warning("⚠️ No matches in test set - may be too specific or incorrect")
                        
            except Exception as e:
                st.error(f"Error processing SMARTS: {str(e)}")
        
        with col2:
            st.write("### Decision")
            
            # Show previous decision if exists
            if current in st.session_state.decisions:
                prev_decision = st.session_state.decisions[current]
                st.info(f"Current: **{prev_decision}**")
            
            # Decision buttons
            col_a, col_b = st.columns(2)
            with col_a:
                if st.button("✅ APPROVE", type="primary", use_container_width=True):
                    st.session_state.decisions[current] = "APPROVED"
                    if current < total - 1:
                        st.session_state.current_idx += 1
                    st.rerun()
            
            with col_b:
                if st.button("❌ REJECT", type="secondary", use_container_width=True):
                    st.session_state.decisions[current] = "REJECTED"
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
                                      value=current + 1, step=1)
            if st.button("Go"):
                st.session_state.current_idx = jump_to - 1
                st.rerun()
    
    else:
        st.success("🎉 You've reviewed all SMARTS patterns!")
    
    # Export results
    if len(st.session_state.decisions) > 0:
        st.write("---")
        st.subheader("Export Results")
        
        # Create results dataframe
        results_df = df.copy()
        results_df['Decision'] = results_df.index.map(
            lambda x: st.session_state.decisions.get(x, "NOT_REVIEWED")
        )
        
        approved = results_df[results_df['Decision'] == 'APPROVED']
        rejected = results_df[results_df['Decision'] == 'REJECTED']
        
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Approved", len(approved))
        with col2:
            st.metric("Rejected", len(rejected))
        with col3:
            st.metric("Pending", total - len(st.session_state.decisions))
        
        # Download buttons
        st.download_button(
            "📥 Download Approved SMARTS",
            approved.to_csv(index=False),
            "approved_smarts.csv",
            "text/csv"
        )
        
        st.download_button(
            "📥 Download Full Results",
            results_df.to_csv(index=False),
            "smarts_review_results.csv",
            "text/csv"
        )

else:
    st.info("👆 Upload a CSV file with your SMARTS patterns to begin review")
    st.write("**Expected format:**")
    st.code("""SMARTS,Description
[#6]-[#8],Carbon-oxygen single bond
c1ccccc1,Benzene ring
[NX3;H2,H1;!$(NC=O)],Primary or secondary amine""")
