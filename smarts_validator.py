import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw, PandasTools
import io
import requests
import base64
import time
import json
from datetime import datetime

# Page config
st.set_page_config(
    page_title="SMARTS Toolkit Pro", 
    layout="wide", 
    initial_sidebar_state="collapsed",
    page_icon="🧬"
)

# ============================================================================
# PROFESSIONAL DESIGN SYSTEM
# ============================================================================

st.markdown("""
<style>
    /* Global Styles */
    .block-container {
        padding-top: 1rem;
        padding-bottom: 0rem;
        padding-left: 2rem;
        padding-right: 2rem;
    }
    
    /* Professional Status Badges */
    .status-badge {
        display: inline-block;
        padding: 0.25rem 0.75rem;
        border-radius: 1rem;
        font-size: 0.75rem;
        font-weight: 600;
        text-transform: uppercase;
        letter-spacing: 0.5px;
    }
    
    .status-pass {
        background: #D1FAE5;
        color: #065F46;
        border: 1px solid #10B981;
    }
    
    .status-review {
        background: #FEF3C7;
        color: #92400E;
        border: 1px solid #F59E0B;
    }
    
    .status-concern {
        background: #FED7AA;
        color: #9A3412;
        border: 1px solid #F97316;
    }
    
    .status-block {
        background: #FEE2E2;
        color: #991B1B;
        border: 1px solid #EF4444;
    }
    
    /* Metric Cards */
    .metric-card {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        padding: 1.5rem;
        border-radius: 0.75rem;
        color: white;
        box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        margin-bottom: 1rem;
    }
    
    .metric-card-value {
        font-size: 2rem;
        font-weight: 700;
        margin: 0.5rem 0;
    }
    
    .metric-card-label {
        font-size: 0.875rem;
        opacity: 0.9;
        text-transform: uppercase;
        letter-spacing: 1px;
    }
    
    /* Dashboard Cards */
    .dashboard-card {
        background: white;
        border: 1px solid #E5E7EB;
        border-radius: 0.5rem;
        padding: 1.5rem;
        margin-bottom: 1rem;
        box-shadow: 0 1px 3px rgba(0,0,0,0.05);
        transition: box-shadow 0.2s;
    }
    
    .dashboard-card:hover {
        box-shadow: 0 4px 6px rgba(0,0,0,0.1);
    }
    
    .dashboard-card-title {
        font-size: 1.125rem;
        font-weight: 600;
        margin-bottom: 0.5rem;
        color: #1F2937;
    }
    
    .dashboard-card-subtitle {
        font-size: 0.875rem;
        color: #6B7280;
        margin-bottom: 1rem;
    }
    
    /* Progress Bar */
    .custom-progress {
        background: #E5E7EB;
        border-radius: 1rem;
        height: 0.5rem;
        overflow: hidden;
    }
    
    .custom-progress-bar {
        background: linear-gradient(90deg, #667eea 0%, #764ba2 100%);
        height: 100%;
        transition: width 0.3s ease;
    }
    
    /* Compact Metrics */
    [data-testid="stMetric"] {
        background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
        padding: 0.75rem;
        border-radius: 0.5rem;
        border: 1px solid #E5E7EB;
    }
    
    [data-testid="stMetricLabel"] {
        font-size: 0.75rem;
        font-weight: 600;
        text-transform: uppercase;
        letter-spacing: 0.5px;
    }
    
    [data-testid="stMetricValue"] {
        font-size: 1.5rem;
        font-weight: 700;
    }
    
    /* Typography */
    h1 {
        font-size: 2rem;
        font-weight: 700;
        margin-bottom: 0.5rem;
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
    }
    
    h2 {
        font-size: 1.5rem;
        font-weight: 600;
        margin-bottom: 0.5rem;
        color: #1F2937;
    }
    
    h3 {
        font-size: 1.25rem;
        font-weight: 600;
        margin-bottom: 0.5rem;
        color: #374151;
    }
    
    /* Toast Notifications */
    .toast-success {
        background: #D1FAE5;
        color: #065F46;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #10B981;
        margin-bottom: 1rem;
    }
    
    .toast-warning {
        background: #FEF3C7;
        color: #92400E;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #F59E0B;
        margin-bottom: 1rem;
    }
    
    /* Simple Charts */
    .chart-bar {
        background: linear-gradient(90deg, #667eea 0%, #764ba2 100%);
        height: 2rem;
        border-radius: 0.25rem;
        margin-bottom: 0.5rem;
        display: flex;
        align-items: center;
        padding: 0 1rem;
        color: white;
        font-weight: 600;
        font-size: 0.875rem;
    }
</style>
""", unsafe_allow_html=True)

# ============================================================================
# REUSABLE COMPONENTS
# ============================================================================

def status_badge(status):
    """Generate professional status badge"""
    badge_map = {
        "OK": ("PASS", "status-pass"),
        "Yellow": ("REVIEW", "status-review"),
        "Amber": ("CONCERN", "status-concern"),
        "Red": ("BLOCK", "status-block"),
        "FLAGGED": ("FLAGGED", "status-block"),
        "CHECKED": ("CHECKED", "status-pass")
    }
    
    label, css_class = badge_map.get(status, (status, "status-pass"))
    return f'<span class="status-badge {css_class}">{label}</span>'

def metric_card(label, value):
    """Generate gradient metric card"""
    return f"""
    <div class="metric-card">
        <div class="metric-card-label">{label}</div>
        <div class="metric-card-value">{value}</div>
    </div>
    """

def progress_bar(current, total):
    """Generate custom progress bar"""
    percentage = (current / total * 100) if total > 0 else 0
    return f"""
    <div class="custom-progress">
        <div class="custom-progress-bar" style="width: {percentage}%"></div>
    </div>
    <p style="text-align: center; font-size: 0.875rem; color: #6B7280; margin-top: 0.5rem;">
        {current} / {total} patterns reviewed ({percentage:.1f}%)
    </p>
    """

def simple_bar_chart(data_dict, max_width=300):
    """Generate simple HTML bar chart"""
    if not data_dict:
        return ""
    
    max_value = max(data_dict.values()) if data_dict.values() else 1
    
    html = ""
    for label, value in data_dict.items():
        width = (value / max_value * max_width) if max_value > 0 else 0
        html += f"""
        <div style="margin-bottom: 0.5rem;">
            <div style="font-size: 0.875rem; color: #6B7280; margin-bottom: 0.25rem;">
                {status_badge(label)} {value}
            </div>
            <div class="chart-bar" style="width: {width}px;">
                {value}
            </div>
        </div>
        """
    return html

# ============================================================================
# SESSION MANAGEMENT
# ============================================================================

def save_session():
    """Save current session to file"""
    session_data = {
        'timestamp': datetime.now().isoformat(),
        'decisions': st.session_state.decisions,
        'current_idx': st.session_state.current_idx,
        'mode': st.session_state.mode,
        'uploaded_file_name': st.session_state.uploaded_file_name
    }
    
    session_file = f"session_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    return json.dumps(session_data, indent=2), session_file

def load_session(uploaded_file):
    """Load session from file"""
    try:
        session_data = json.load(uploaded_file)
        st.session_state.decisions = session_data.get('decisions', {})
        st.session_state.decisions = {int(k): v for k, v in st.session_state.decisions.items()}
        st.session_state.current_idx = session_data.get('current_idx', 0)
        return True
    except Exception as e:
        st.error(f"Failed to load session: {str(e)}")
        return False

# ============================================================================
# ANALYTICS FUNCTIONS
# ============================================================================

def calculate_session_stats():
    """Calculate comprehensive session statistics"""
    if not st.session_state.decisions:
        return None
    
    total_reviewed = len(st.session_state.decisions)
    decisions_count = {}
    for decision in st.session_state.decisions.values():
        decisions_count[decision] = decisions_count.get(decision, 0) + 1
    
    stats = {
        'total_reviewed': total_reviewed,
        'pass_rate': decisions_count.get('OK', 0) / total_reviewed * 100 if total_reviewed > 0 else 0,
        'block_rate': decisions_count.get('Red', 0) / total_reviewed * 100 if total_reviewed > 0 else 0,
        'decisions_count': decisions_count
    }
    
    return stats

# ============================================================================
# INITIALIZE SESSION STATE
# ============================================================================

if 'mode' not in st.session_state:
    st.session_state.mode = "Dashboard"
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
if 'search_query' not in st.session_state:
    st.session_state.search_query = ""
if 'active_filters' not in st.session_state:
    st.session_state.active_filters = []
if 'session_start_time' not in st.session_state:
    st.session_state.session_start_time = datetime.now()
if 'session_name' not in st.session_state:
    st.session_state.session_name = f"Session_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
if 'quick_test_results' not in st.session_state:
    st.session_state.quick_test_results = None
if 'quick_page_idx' not in st.session_state:
    st.session_state.quick_page_idx = 0

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def get_smartsplus_image(smarts_pattern, api_key, use_cache=True):
    """Get visualization from SMARTS.plus API"""
    cache_key = f"{smarts_pattern}_{api_key[:8]}_svg"
    if use_cache and cache_key in st.session_state.viz_cache:
        return st.session_state.viz_cache[cache_key]
    
    try:
        payload = {
            "query": {
                "smarts": smarts_pattern,
                "parameters": {
                    "file_format": "svg",
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
        if st.session_state.get('debug_mode', False):
            st.error(f"SMARTS.plus error: {str(e)}")
        return None

def filter_patterns(df, search_query, active_filters):
    """Filter patterns based on search and filters"""
    if df is None:
        return None
    
    filtered_df = df.copy()
    
    if search_query:
        mask = (
            filtered_df['SMARTS'].str.contains(search_query, case=False, na=False) |
            (filtered_df['Description'].str.contains(search_query, case=False, na=False) if 'Description' in filtered_df.columns else False)
        )
        filtered_df = filtered_df[mask]
    
    if "Not Reviewed" in active_filters:
        reviewed_indices = set(st.session_state.decisions.keys())
        filtered_df = filtered_df[~filtered_df.index.isin(reviewed_indices)]
    
    if "Flagged Only" in active_filters:
        flagged_indices = [idx for idx, dec in st.session_state.decisions.items() if dec in ["Red", "Amber"]]
        filtered_df = filtered_df[filtered_df.index.isin(flagged_indices)]
    
    return filtered_df

# ============================================================================
# HEADER & NAVIGATION
# ============================================================================

st.title("🧬 SMARTS Toolkit Pro")

col_nav1, col_nav2, col_nav3 = st.columns([1, 3, 1])

with col_nav1:
    st.markdown(f"**Session:** {st.session_state.session_name}")

with col_nav2:
    mode = st.radio(
        "Navigation:",
        ["Dashboard", "Pattern Curator", "Batch Validator", "Quick Tester", "Production Filter", "Analytics"],
        horizontal=True,
        label_visibility="collapsed"
    )
    st.session_state.mode = mode

with col_nav3:
    duration = datetime.now() - st.session_state.session_start_time
    st.markdown(f"**Time:** {str(duration).split('.')[0]}")

st.divider()

# Settings expander
with st.expander("⚙️ Settings", expanded=False):
    col_s1, col_s2, col_s3 = st.columns(3)
    
    with col_s1:
        st.subheader("API Settings")
        api_key_input = st.text_input(
            "SMARTS.plus API Key",
            type="password",
            value=st.session_state.api_key,
            help="Get free key at smarts.plus/sign_up"
        )
        if api_key_input != st.session_state.api_key:
            st.session_state.api_key = api_key_input
            st.session_state.viz_cache = {}
        
        use_smartsplus = st.checkbox(
            "Use SMARTS.plus visualization",
            value=bool(st.session_state.api_key),
            disabled=not bool(st.session_state.api_key)
        )
    
    with col_s2:
        st.subheader("Session Management")
        
        session_name_input = st.text_input(
            "Session Name",
            value=st.session_state.session_name
        )
        st.session_state.session_name = session_name_input
        
        col_save, col_load = st.columns(2)
        with col_save:
            if st.button("💾 Save", use_container_width=True):
                session_json, filename = save_session()
                st.download_button(
                    "📥 Download",
                    session_json,
                    filename,
                    "application/json",
                    use_container_width=True
                )
        
        with col_load:
            uploaded_session = st.file_uploader(
                "Load Session",
                type=['json'],
                key='session_upload',
                label_visibility="collapsed"
            )
            if uploaded_session:
                if load_session(uploaded_session):
                    st.success("Loaded!")
                    st.rerun()
    
    with col_s3:
        st.subheader("Info")
        st.session_state.debug_mode = st.checkbox("Debug mode", value=False)
        st.metric("Reviewed", len(st.session_state.decisions))
        st.metric("Cache", len(st.session_state.viz_cache))

st.divider()

# ============================================================================
# MODE: DASHBOARD
# ============================================================================

if mode == "Dashboard":
    st.header("📊 Dashboard")
    
    stats = calculate_session_stats()
    
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        total = stats['total_reviewed'] if stats else 0
        st.markdown(metric_card("Reviewed", total), unsafe_allow_html=True)
    
    with col2:
        pass_rate = stats['pass_rate'] if stats else 0
        st.markdown(metric_card("Pass Rate", f"{pass_rate:.1f}%"), unsafe_allow_html=True)
    
    with col3:
        block_rate = stats['block_rate'] if stats else 0
        st.markdown(metric_card("Block Rate", f"{block_rate:.1f}%"), unsafe_allow_html=True)
    
    with col4:
        total_patterns = len(st.session_state.smarts_data) if st.session_state.smarts_data is not None else 0
        st.markdown(metric_card("Total", total_patterns), unsafe_allow_html=True)
    
    st.write("")
    
    if stats:
        st.subheader("Decision Distribution")
        chart_html = simple_bar_chart(stats['decisions_count'])
        st.markdown(chart_html, unsafe_allow_html=True)
    
    st.write("")
    st.subheader("Quick Actions")
    
    col_a1, col_a2, col_a3 = st.columns(3)
    
    with col_a1:
        st.markdown("""
        <div class="dashboard-card">
            <div class="dashboard-card-title">🔬 Pattern Curator</div>
            <div class="dashboard-card-subtitle">Review and classify patterns</div>
        </div>
        """, unsafe_allow_html=True)
        if st.button("Start Curation", use_container_width=True, type="primary"):
            st.session_state.mode = "Pattern Curator"
            st.rerun()
    
    with col_a2:
        st.markdown("""
        <div class="dashboard-card">
            <div class="dashboard-card-title">⚡ Quick Tester</div>
            <div class="dashboard-card-subtitle">Test patterns instantly</div>
        </div>
        """, unsafe_allow_html=True)
        if st.button("Quick Test", use_container_width=True):
            st.session_state.mode = "Quick Tester"
            st.rerun()
    
    with col_a3:
        st.markdown("""
        <div class="dashboard-card">
            <div class="dashboard-card-title">📈 Analytics</div>
            <div class="dashboard-card-subtitle">View insights</div>
        </div>
        """, unsafe_allow_html=True)
        if st.button("Analytics", use_container_width=True):
            st.session_state.mode = "Analytics"
            st.rerun()

# ============================================================================
# MODE: PATTERN CURATOR
# ============================================================================

elif mode == "Pattern Curator":
    st.header("🔬 Pattern Curator")
    
    col_search, col_filter = st.columns([3, 2])
    
    with col_search:
        search_query = st.text_input(
            "🔍 Search",
            placeholder="Search patterns...",
            key="curator_search",
            label_visibility="collapsed"
        )
        st.session_state.search_query = search_query
    
    with col_filter:
        active_filters = st.multiselect(
            "Filters",
            ["Not Reviewed", "Flagged Only"],
            key="curator_filters",
            label_visibility="collapsed"
        )
        st.session_state.active_filters = active_filters
    
    st.write("")
    
    uploaded_file = st.file_uploader("📁 Upload SMARTS CSV", type=['csv'], label_visibility="collapsed")
    
    if uploaded_file:
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
                st.success(f"✅ Loaded {len(df)} patterns")
            except Exception as e:
                st.error(f"Error: {str(e)}")
                st.stop()
        
        df = st.session_state.smarts_data
        filtered_df = filter_patterns(df, search_query, active_filters)
        
        if filtered_df is None or len(filtered_df) == 0:
            st.warning("No patterns match")
            st.stop()
        
        filtered_indices = filtered_df.index.tolist()
        
        if st.session_state.current_idx in filtered_indices:
            current_filtered_pos = filtered_indices.index(st.session_state.current_idx)
        else:
            current_filtered_pos = 0
            st.session_state.current_idx = filtered_indices[0]
        
        total = len(filtered_df)
        current = current_filtered_pos
        actual_idx = filtered_indices[current]
        
        progress = len(st.session_state.decisions) / len(df) if len(df) > 0 else 0
        st.markdown(progress_bar(len(st.session_state.decisions), len(df)), unsafe_allow_html=True)
        
        st.write("")
        
        col1, col2, col3, col4, col5 = st.columns(5)
        with col1:
            st.metric("Reviewed", f"{len(st.session_state.decisions)}/{len(df)}")
        with col2:
            pass_count = sum(1 for d in st.session_state.decisions.values() if d == "OK")
            st.metric("Pass", pass_count)
        with col3:
            yellow_count = sum(1 for d in st.session_state.decisions.values() if d == "Yellow")
            st.metric("Review", yellow_count)
        with col4:
            amber_count = sum(1 for d in st.session_state.decisions.values() if d == "Amber")
            st.metric("Concern", amber_count)
        with col5:
            red_count = sum(1 for d in st.session_state.decisions.values() if d == "Red")
            st.metric("Block", red_count)
        
        st.write("")
        
        if current < total:
            col_main, col_side = st.columns([2.5, 1])
            
            with col_main:
                smarts_pattern = df.iloc[actual_idx]['SMARTS']
                
                st.markdown(f"""
                <div class="dashboard-card">
                    <div class="dashboard-card-title">Pattern {current + 1} of {total}</div>
                </div>
                """, unsafe_allow_html=True)
                
                st.code(smarts_pattern, language='text')
                
                if 'Description' in df.columns and pd.notna(df.iloc[actual_idx]['Description']):
                    st.caption(df.iloc[actual_idx]['Description'])
                
                try:
                    pattern = Chem.MolFromSmarts(smarts_pattern)
                    if pattern is None:
                        st.error("⚠️ Invalid SMARTS")
                    else:
                        image_displayed = False
                        if use_smartsplus and st.session_state.api_key:
                            with st.spinner("Loading..."):
                                smartsplus_svg = get_smartsplus_image(smarts_pattern, st.session_state.api_key)
                                if smartsplus_svg:
                                    import re
                                    svg_fixed = re.sub(r'\s*width="[^"]*"', '', smartsplus_svg)
                                    svg_fixed = re.sub(r'\s*height="[^"]*"', '', svg_fixed)
                                    st.markdown(f'<div style="width: 800px; max-width: 100%;">{svg_fixed}</div>', unsafe_allow_html=True)
                                    st.caption("🎨 SMARTS.plus")
                                    image_displayed = True
                        
                        if not image_displayed:
                            img = Draw.MolToImage(pattern, size=(400, 320))
                            st.image(img, width=400)
                except Exception as e:
                    st.error(f"Error: {str(e)}")
            
            with col_side:
                if actual_idx in st.session_state.decisions:
                    current_decision = st.session_state.decisions[actual_idx]
                    st.markdown(f"**Current:** {status_badge(current_decision)}", unsafe_allow_html=True)
                    st.write("")
                
                st.markdown("**Classify:**")
                
                if st.button("✓ PASS", use_container_width=True, key="pass", type="primary"):
                    st.session_state.decisions[actual_idx] = "OK"
                    if current < total - 1:
                        st.session_state.current_idx = filtered_indices[current + 1]
                    st.rerun()
                
                if st.button("👁️ REVIEW", use_container_width=True, key="review"):
                    st.session_state.decisions[actual_idx] = "Yellow"
                    if current < total - 1:
                        st.session_state.current_idx = filtered_indices[current + 1]
                    st.rerun()
                
                if st.button("⚠️ CONCERN", use_container_width=True, key="concern"):
                    st.session_state.decisions[actual_idx] = "Amber"
                    if current < total - 1:
                        st.session_state.current_idx = filtered_indices[current + 1]
                    st.rerun()
                
                if st.button("🛑 BLOCK", use_container_width=True, key="block"):
                    st.session_state.decisions[actual_idx] = "Red"
                    if current < total - 1:
                        st.session_state.current_idx = filtered_indices[current + 1]
                    st.rerun()
                
                st.divider()
                
                col_n1, col_n2 = st.columns(2)
                with col_n1:
                    if st.button("⬅️", disabled=(current == 0), use_container_width=True):
                        st.session_state.current_idx = filtered_indices[current - 1]
                        st.rerun()
                with col_n2:
                    if st.button("➡️", disabled=(current >= total - 1), use_container_width=True):
                        st.session_state.current_idx = filtered_indices[current + 1]
                        st.rerun()
                
                jump_to = st.number_input("Jump to", 1, total, current + 1, key='jump', label_visibility="collapsed")
                if st.button("🎯 Jump", use_container_width=True):
                    st.session_state.current_idx = filtered_indices[jump_to - 1]
                    st.rerun()
        
        else:
            st.success("🎉 Complete!")
        
        if len(st.session_state.decisions) > 0:
            st.write("")
            st.divider()
            st.subheader("📥 Export")
            
            results_df = df.copy()
            results_df['Decision'] = results_df.index.map(lambda x: st.session_state.decisions.get(x, "NOT_REVIEWED"))
            
            col_e1, col_e2, col_e3, col_e4 = st.columns(4)
            
            with col_e1:
                st.download_button("📥 All", results_df.to_csv(index=False), "all.csv", "text/csv", use_container_width=True)
            
            with col_e2:
                concerning = results_df[results_df['Decision'].isin(['Amber', 'Red'])]
                if len(concerning) > 0:
                    st.download_button("⚠️ Concern", concerning.to_csv(index=False), "concern.csv", "text/csv", use_container_width=True)
            
            with col_e3:
                ok_df = results_df[results_df['Decision'] == 'OK']
                if len(ok_df) > 0:
                    st.download_button("✅ Pass", ok_df.to_csv(index=False), "pass.csv", "text/csv", use_container_width=True)
            
            with col_e4:
                session_json, filename = save_session()
                st.download_button("💾 Session", session_json, filename, "application/json", use_container_width=True)
    
    else:
        st.info("👆 Upload CSV to begin")

# ============================================================================
# MODE: BATCH VALIDATOR & QUICK TESTER
# ============================================================================

elif mode == "Batch Validator":
    st.header("🔬 Batch Validator")
    
    validator_mode = st.radio("Type:", ["Batch Validation", "Quick Test"], horizontal=True, key="val_submode")
    st.write("")
    
    if validator_mode == "Batch Validation":
        col_u1, col_u2 = st.columns(2)
        with col_u1:
            smarts_file = st.file_uploader("📁 SMARTS", type=['csv'], key='val_smarts', label_visibility="collapsed")
        with col_u2:
            molecules_file = st.file_uploader("📁 Molecules", type=['csv'], key='val_mols', label_visibility="collapsed")
        
        if smarts_file and st.session_state.smarts_data is None:
            try:
                df = pd.read_csv(smarts_file)
                if 'SMARTS' not in df.columns:
                    df.columns = ['SMARTS'] + list(df.columns[1:])
                st.session_state.smarts_data = df
                st.session_state.current_idx = 0
                st.session_state.decisions = {}
                st.success(f"✅ Loaded {len(df)} patterns")
            except Exception as e:
                st.error(f"Error: {str(e)}")
        
        if molecules_file and st.session_state.test_molecules is None:
            try:
                mol_df = pd.read_csv(molecules_file)
                smiles_col = 'SMILES' if 'SMILES' in mol_df.columns else mol_df.columns[0]
                mol_df['Mol'] = mol_df[smiles_col].apply(lambda x: Chem.MolFromSmiles(x) if pd.notna(x) else None)
                mol_df = mol_df[mol_df['Mol'].notna()]
                st.session_state.test_molecules = mol_df
                st.success(f"✅ Loaded {len(mol_df)} molecules")
            except Exception as e:
                st.error(f"Error: {str(e)}")
        
        if st.session_state.smarts_data is not None and st.session_state.test_molecules is not None:
            df = st.session_state.smarts_data
            mol_df = st.session_state.test_molecules
            total = len(df)
            current = st.session_state.current_idx
            
            st.markdown(progress_bar(current, total), unsafe_allow_html=True)
            st.write("")
            
            col_m1, col_m2 = st.columns(2)
            with col_m1:
                flagged = sum(1 for d in st.session_state.decisions.values() if d == "FLAGGED")
                st.metric("Flagged", flagged)
            with col_m2:
                checked = sum(1 for d in st.session_state.decisions.values() if d == "CHECKED")
                st.metric("Checked", checked)
            
            st.write("")
            
            if current < total:
                smarts_pattern = df.iloc[current]['SMARTS']
                
                col_main, col_side = st.columns([2.5, 1])
                
                with col_main:
                    st.markdown(f'<div class="dashboard-card"><div class="dashboard-card-title">Pattern {current + 1} of {total}</div></div>', unsafe_allow_html=True)
                    st.code(smarts_pattern, language='text')
                    
                    if 'Description' in df.columns and pd.notna(df.iloc[current]['Description']):
                        st.caption(df.iloc[current]['Description'])
                    
                    try:
                        pattern = Chem.MolFromSmarts(smarts_pattern)
                        if pattern:
                            matches = []
                            for idx, row in mol_df.iterrows():
                                mol = row.get('Mol')
                                if mol and mol.HasSubstructMatch(pattern):
                                    matches.append((idx, mol))
                            
                            match_rate = len(matches)/len(mol_df)*100
                            st.markdown(f"**{len(matches)}/{len(mol_df)} matches ({match_rate:.1f}%)**")
                            
                            if match_rate > 50:
                                st.warning("⚠️ High match rate")
                            elif match_rate == 0:
                                st.info("ℹ️ No matches")
                            
                            if len(matches) > 0:
                                st.write("**Samples:**")
                                cols = st.columns(3)
                                for i, (idx, mol) in enumerate(matches[:6]):
                                    with cols[i % 3]:
                                        img = Draw.MolToImage(mol, size=(150, 150), highlightAtoms=mol.GetSubstructMatch(pattern))
                                        st.image(img, width=150)
                                        if 'Name' in mol_df.columns:
                                            st.caption(mol_df.iloc[idx]['Name'])
                    except Exception as e:
                        st.error(f"Error: {str(e)}")
                
                with col_side:
                    if current in st.session_state.decisions:
                        st.markdown(status_badge(st.session_state.decisions[current]), unsafe_allow_html=True)
                    
                    st.write("")
                    
                    if st.button("🚩 FLAG", use_container_width=True, type="primary"):
                        st.session_state.decisions[current] = "FLAGGED"
                        st.rerun()
                    
                    if st.button("✓ OK", use_container_width=True):
                        st.session_state.decisions[current] = "CHECKED"
                        if current < total - 1:
                            st.session_state.current_idx += 1
                        st.rerun()
                    
                    st.divider()
                    
                    col_n1, col_n2 = st.columns(2)
                    with col_n1:
                        if st.button("⬅️", disabled=(current == 0), use_container_width=True):
                            st.session_state.current_idx -= 1
                            st.rerun()
                    with col_n2:
                        if st.button("➡️", disabled=(current >= total - 1), use_container_width=True):
                            st.session_state.current_idx += 1
                            st.rerun()
            else:
                st.success("🎉 Complete!")
            
            if len(st.session_state.decisions) > 0:
                st.write("")
                st.divider()
                st.subheader("📥 Export")
                
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
            st.info("👆 Upload both files to begin")

elif mode == "Quick Tester":
    st.header("⚡ Quick Tester")
    
    col_input, col_output = st.columns([1, 1])
    
    with col_input:
        st.subheader("Input")
        
        tab1, tab2 = st.tabs(["Custom", "Library"])
        
        with tab1:
            quick_smarts = st.text_area("SMARTS:", height=100, placeholder="c1ccccc1", key="quick_smarts_custom")
        
        with tab2:
            st.write("**Common Patterns:**")
            pattern_library = {
                "Benzene": "c1ccccc1",
                "Primary Amine": "[NX3;H2;!$(NC=O)]",
                "Carboxylic Acid": "C(=O)[OH]",
                "Aromatic": "a",
                "Ester": "C(=O)O",
                "Amide": "C(=O)N",
                "Ketone": "[CX3](=O)[#6]",
                "Aldehyde": "[CX3H1](=O)[#6]",
                "Ether": "[OD2]([#6])[#6]"
            }
            
            selected = st.selectbox("Select:", list(pattern_library.keys()))
            if st.button("Use Pattern", use_container_width=True):
                quick_smarts = pattern_library[selected]
                st.session_state.quick_smarts_custom = quick_smarts
                st.rerun()
        
        if quick_smarts:
            pattern = Chem.MolFromSmarts(quick_smarts)
            if pattern:
                st.success("✓ Valid")
            else:
                st.error("✗ Invalid")
        
        st.divider()
        
        mol_source = st.radio("Test:", ["Single SMILES", "File"], horizontal=True)
        
        if mol_source == "Single SMILES":
            quick_smiles = st.text_input("SMILES:", placeholder="CCO")
            
            if quick_smarts and quick_smiles:
                if st.button("🚀 Test", type="primary", use_container_width=True):
                    try:
                        pattern = Chem.MolFromSmarts(quick_smarts)
                        mol = Chem.MolFromSmiles(quick_smiles)
                        
                        if pattern and mol:
                            is_match = mol.HasSubstructMatch(pattern)
                            st.session_state.quick_test_results = {
                                'is_match': is_match,
                                'mol': mol,
                                'pattern': pattern
                            }
                        else:
                            st.error("Invalid input")
                    except Exception as e:
                        st.error(f"Error: {str(e)}")
        else:
            quick_mol_file = st.file_uploader("📁 Molecules", type=['csv'], key='quick_mols', label_visibility="collapsed")
            
            if quick_mol_file and quick_smarts:
                if st.button("🚀 Run", type="primary", use_container_width=True):
                    try:
                        mol_df = pd.read_csv(quick_mol_file)
                        smiles_col = 'SMILES' if 'SMILES' in mol_df.columns else mol_df.columns[0]
                        mol_df['Mol'] = mol_df[smiles_col].apply(lambda x: Chem.MolFromSmiles(x) if pd.notna(x) else None)
                        mol_df = mol_df[mol_df['Mol'].notna()]
                        
                        pattern = Chem.MolFromSmarts(quick_smarts)
                        if pattern:
                            matches = []
                            for idx, row in mol_df.iterrows():
                                mol = row.get('Mol')
                                if mol and mol.HasSubstructMatch(pattern):
                                    matches.append((idx, mol, row))
                            
                            st.session_state.quick_test_results = {
                                'pattern': pattern,
                                'matches': matches,
                                'total': len(mol_df),
                                'mol_df': mol_df
                            }
                            st.session_state.quick_page_idx = 0
                    except Exception as e:
                        st.error(f"Error: {str(e)}")
    
    with col_output:
        st.subheader("Results")
        
        if st.session_state.quick_test_results:
            results = st.session_state.quick_test_results
            
            if 'is_match' in results:
                is_match = results['is_match']
                mol = results['mol']
                pattern = results['pattern']
                
                if is_match:
                    st.markdown('<div class="toast-success"><strong>✅ MATCH</strong></div>', unsafe_allow_html=True)
                    match_atoms = mol.GetSubstructMatch(pattern)
                    img = Draw.MolToImage(mol, size=(400, 400), highlightAtoms=match_atoms)
                else:
                    st.markdown('<div class="toast-warning"><strong>❌ NO MATCH</strong></div>', unsafe_allow_html=True)
                    img = Draw.MolToImage(mol, size=(400, 400))
                
                st.image(img, use_column_width=True)
            
            elif 'matches' in results:
                matches = results['matches']
                total = results['total']
                
                col_r1, col_r2, col_r3 = st.columns(3)
                with col_r1:
                    st.metric("Total", total)
                with col_r2:
                    st.metric("Matches", len(matches))
                with col_r3:
                    rate = len(matches)/total*100 if total > 0 else 0
                    st.metric("Rate", f"{rate:.1f}%")
                
                st.write("")
                
                if len(matches) > 0:
                    per_page = 6
                    total_pages = (len(matches) + per_page - 1) // per_page
                    current_page = st.session_state.quick_page_idx
                    
                    start = current_page * per_page
                    end = min(start + per_page, len(matches))
                    page_matches = matches[start:end]
                    
                    st.write(f"**Showing {start + 1}-{end} of {len(matches)}**")
                    
                    cols = st.columns(3)
                    for i, (idx, mol, row) in enumerate(page_matches):
                        with cols[i % 3]:
                            match_atoms = mol.GetSubstructMatch(results['pattern'])
                            img = Draw.MolToImage(mol, size=(150, 150), highlightAtoms=match_atoms)
                            st.image(img, width=150)
                            if 'Name' in row:
                                st.caption(row['Name'])
                    
                    if total_pages > 1:
                        st.write("")
                        col_p1, col_p2, col_p3 = st.columns([1, 2, 1])
                        with col_p1:
                            if st.button("⬅️", disabled=(current_page == 0), use_container_width=True):
                                st.session_state.quick_page_idx -= 1
                                st.rerun()
                        with col_p2:
                            st.write(f"Page {current_page + 1}/{total_pages}")
                        with col_p3:
                            if st.button("➡️", disabled=(current_page >= total_pages - 1), use_container_width=True):
                                st.session_state.quick_page_idx += 1
                                st.rerun()
                    
                    st.write("")
                    match_indices = [idx for idx, _, _ in matches]
                    export_df = results['mol_df'].iloc[match_indices].copy()
                    cols = [c for c in export_df.columns if c not in ['Mol']]
                    if 'SMILES' not in cols:
                        export_df['SMILES'] = export_df['Mol'].apply(lambda m: Chem.MolToSmiles(m) if m else '')
                        cols = ['SMILES'] + cols
                    st.download_button("📥 Download", export_df[cols].to_csv(index=False), "matches.csv", "text/csv", use_container_width=True, type="primary")
                else:
                    st.info("No matches")
        else:
            st.info("👈 Configure and run test")

# ============================================================================
# MODE: PRODUCTION FILTER
# ============================================================================

elif mode == "Production Filter":
    st.header("🚀 Production Filter")
    
    col_u1, col_u2 = st.columns(2)
    with col_u1:
        filter_file = st.file_uploader("📁 Approved SMARTS", type=['csv'], key='filter_smarts')
    with col_u2:
        gen_ai_file = st.file_uploader("📁 Gen AI Molecules", type=['csv'], key='gen_ai_mols')
    
    filter_severity = st.multiselect("Apply:", ["Red", "Amber", "Yellow", "OK"], default=["Red", "Amber"])
    
    filter_patterns = None
    if filter_file:
        try:
            filter_df = pd.read_csv(filter_file)
            if 'Decision' in filter_df.columns:
                filter_patterns = filter_df[filter_df['Decision'].isin(filter_severity)]
                st.success(f"✅ {len(filter_patterns)} patterns")
            else:
                st.warning("No 'Decision' column")
                filter_patterns = filter_df
        except Exception as e:
            st.error(f"Error: {str(e)}")
    
    gen_ai_mols = None
    if gen_ai_file:
        try:
            gen_ai_mols = pd.read_csv(gen_ai_file)
            smiles_col = 'SMILES' if 'SMILES' in gen_ai_mols.columns else gen_ai_mols.columns[0]
            gen_ai_mols['Mol'] = gen_ai_mols[smiles_col].apply(lambda x: Chem.MolFromSmiles(x) if pd.notna(x) else None)
            gen_ai_mols = gen_ai_mols[gen_ai_mols['Mol'].notna()]
            st.success(f"✅ {len(gen_ai_mols)} molecules")
        except Exception as e:
            st.error(f"Error: {str(e)}")
    
    if filter_patterns is not None and gen_ai_mols is not None and len(filter_patterns) > 0:
        if st.button("🚀 Run Filter", type="primary", use_container_width=True):
            with st.spinner("Filtering..."):
                compiled = []
                for idx, row in filter_patterns.iterrows():
                    p = Chem.MolFromSmarts(row['SMARTS'])
                    if p:
                        compiled.append({
                            'pattern': p,
                            'smarts': row['SMARTS'],
                            'description': row.get('Description', ''),
                            'decision': row.get('Decision', 'Unknown')
                        })
                
                passed = []
                failed = []
                rejection_reasons = []
                
                progress_placeholder = st.empty()
                
                for idx, row in gen_ai_mols.iterrows():
                    mol = row.get('Mol')
                    if mol is None:
                        continue
                    
                    progress = (idx + 1) / len(gen_ai_mols)
                    progress_placeholder.progress(progress, text=f"Processing {idx + 1}/{len(gen_ai_mols)}")
                    
                    caught_by = []
                    for p in compiled:
                        if mol.HasSubstructMatch(p['pattern']):
                            caught_by.append(p)
                    
                    if len(caught_by) == 0:
                        passed.append(idx)
                    else:
                        failed.append(idx)
                        rejection_reasons.append({
                            'molecule_idx': idx,
                            'num_violations': len(caught_by),
                            'patterns': [p['smarts'] for p in caught_by],
                            'decisions': [p['decision'] for p in caught_by]
                        })
                
                progress_placeholder.empty()
                
                st.session_state.filter_results = {
                    'passed': passed,
                    'failed': failed,
                    'rejection_reasons': rejection_reasons,
                    'total': len(gen_ai_mols),
                    'patterns_used': len(compiled)
                }
                st.rerun()
    
    if st.session_state.filter_results:
        results = st.session_state.filter_results
        
        st.write("")
        st.divider()
        st.subheader("📊 Results")
        
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.markdown(metric_card("Total", results['total']), unsafe_allow_html=True)
        with col2:
            pass_rate = len(results['passed'])/results['total']*100
            st.markdown(metric_card("Passed", f"{len(results['passed'])} ({pass_rate:.1f}%)"), unsafe_allow_html=True)
        with col3:
            fail_rate = len(results['failed'])/results['total']*100
            st.markdown(metric_card("Failed", f"{len(results['failed'])} ({fail_rate:.1f}%)"), unsafe_allow_html=True)
        with col4:
            st.markdown(metric_card("Filters", results['patterns_used']), unsafe_allow_html=True)
        
        st.write("")
        
        col_e1, col_e2 = st.columns(2)
        with col_e1:
            if len(results['passed']) > 0:
                passed_mols = gen_ai_mols.iloc[results['passed']].copy()
                if 'SMILES' in passed_mols.columns:
                    export_df = passed_mols[['SMILES']].copy()
                else:
                    export_df = pd.DataFrame({'SMILES': passed_mols['Mol'].apply(lambda m: Chem.MolToSmiles(m) if m else '')})
                st.download_button("📥 Passed", export_df.to_csv(index=False), "passed.csv", "text/csv", use_container_width=True, type="primary")
        
        with col_e2:
            if len(results['failed']) > 0:
                rejection_df = pd.DataFrame(results['rejection_reasons'])
                st.download_button("📥 Rejected", rejection_df.to_csv(index=False), "rejected.csv", "text/csv", use_container_width=True)
        
        if len(results['failed']) > 0:
            with st.expander(f"🔍 Sample Rejections (5 of {len(results['failed'])})"):
                for i, reason in enumerate(results['rejection_reasons'][:5]):
                    st.markdown(f'<div class="dashboard-card"><div class="dashboard-card-title">Molecule {reason["molecule_idx"]}</div><div class="dashboard-card-subtitle">Caught by {reason["num_violations"]} patterns</div></div>', unsafe_allow_html=True)
                    for pattern, decision in zip(reason['patterns'], reason['decisions']):
                        col_r1, col_r2 = st.columns([3, 1])
                        with col_r1:
                            st.code(pattern)
                        with col_r2:
                            st.markdown(status_badge(decision), unsafe_allow_html=True)
                    if i < 4:
                        st.divider()
    elif filter_patterns is not None and gen_ai_mols is not None:
        if len(filter_patterns) == 0:
            st.warning("⚠️ No patterns selected")
    else:
        st.info("👆 Upload files and run filter")

# ============================================================================
# MODE: ANALYTICS
# ============================================================================

elif mode == "Analytics":
    st.header("📈 Analytics")
    
    stats = calculate_session_stats()
    
    if stats:
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.markdown(metric_card("Reviewed", stats['total_reviewed']), unsafe_allow_html=True)
        with col2:
            st.markdown(metric_card("Pass Rate", f"{stats['pass_rate']:.1f}%"), unsafe_allow_html=True)
        with col3:
            st.markdown(metric_card("Block Rate", f"{stats['block_rate']:.1f}%"), unsafe_allow_html=True)
        with col4:
            duration = datetime.now() - st.session_state.session_start_time
            minutes = int(duration.total_seconds() / 60)
            st.markdown(metric_card("Duration", f"{minutes} min"), unsafe_allow_html=True)
        
        st.write("")
        
        col_chart1, col_chart2 = st.columns(2)
        
        with col_chart1:
            st.subheader("Decision Distribution")
            chart_html = simple_bar_chart(stats['decisions_count'], max_width=400)
            st.markdown(chart_html, unsafe_allow_html=True)
        
        with col_chart2:
            st.subheader("Breakdown")
            decisions = stats['decisions_count']
            for decision, count in decisions.items():
                pct = count/stats['total_reviewed']*100
                st.markdown(f"{status_badge(decision)} {count} ({pct:.1f}%)", unsafe_allow_html=True)
                st.write("")
        
        st.write("")
        st.divider()
        st.subheader("🎯 Insights")
        
        if stats['block_rate'] > 20:
            st.markdown('<div class="toast-warning"><strong>⚠️ High Block Rate</strong><br>Over 20% blocked. Review criteria.</div>', unsafe_allow_html=True)
        
        if stats['pass_rate'] > 80:
            st.markdown('<div class="toast-success"><strong>✅ High Pass Rate</strong><br>Great! Over 80% passing.</div>', unsafe_allow_html=True)
        
        if stats['total_reviewed'] < 50:
            st.info("💡 Review at least 50 patterns for meaningful statistics")
    else:
        st.info("No data yet. Start reviewing in Pattern Curator!")

# Footer
st.write("")
st.divider()
st.caption("🧬 SMARTS Toolkit Pro | Built for Digital Chemistry")
