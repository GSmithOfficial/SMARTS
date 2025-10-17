import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw, PandasTools
import io
import requests
import base64
import time
import json
from datetime import datetime, timedelta
import plotly.express as px
import plotly.graph_objects as go

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
    
    .metric-card-delta {
        font-size: 0.875rem;
        margin-top: 0.5rem;
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
    
    /* Search Bar */
    .search-container {
        background: white;
        border: 2px solid #E5E7EB;
        border-radius: 0.5rem;
        padding: 0.75rem;
        margin-bottom: 1rem;
        transition: border-color 0.2s;
    }
    
    .search-container:focus-within {
        border-color: #667eea;
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
    
    /* Action Buttons */
    .action-button-primary {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        border: none;
        padding: 0.75rem 1.5rem;
        border-radius: 0.5rem;
        font-weight: 600;
        cursor: pointer;
        transition: transform 0.2s;
    }
    
    .action-button-primary:hover {
        transform: translateY(-2px);
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
    
    /* Buttons */
    .stButton button {
        border-radius: 0.5rem;
        font-weight: 500;
        transition: all 0.2s;
    }
    
    .stButton button:hover {
        transform: translateY(-1px);
        box-shadow: 0 4px 6px rgba(0,0,0,0.1);
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
    
    /* Session Info */
    .session-info {
        background: #F3F4F6;
        padding: 1rem;
        border-radius: 0.5rem;
        font-size: 0.875rem;
        color: #6B7280;
        margin-bottom: 1rem;
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

def metric_card(label, value, delta=None, color="purple"):
    """Generate gradient metric card"""
    delta_html = f'<div class="metric-card-delta">{"↑" if delta and delta > 0 else "↓"} {abs(delta) if delta else ""}</div>' if delta else ""
    
    return f"""
    <div class="metric-card">
        <div class="metric-card-label">{label}</div>
        <div class="metric-card-value">{value}</div>
        {delta_html}
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
        # Convert string keys back to integers
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

def generate_analytics_charts():
    """Generate analytics visualizations"""
    stats = calculate_session_stats()
    if not stats:
        return None, None
    
    # Decision distribution pie chart
    decision_data = stats['decisions_count']
    
    fig_pie = go.Figure(data=[go.Pie(
        labels=list(decision_data.keys()),
        values=list(decision_data.values()),
        marker=dict(colors=['#10B981', '#F59E0B', '#F97316', '#EF4444']),
        hole=0.4
    )])
    fig_pie.update_layout(
        title="Decision Distribution",
        height=300,
        margin=dict(l=20, r=20, t=40, b=20)
    )
    
    # Review velocity (mock data for now - could track real timestamps)
    review_velocity = go.Figure(data=[go.Scatter(
        x=list(range(len(st.session_state.decisions))),
        y=list(range(len(st.session_state.decisions))),
        mode='lines+markers',
        line=dict(color='#667eea', width=2)
    )])
    review_velocity.update_layout(
        title="Review Progress",
        xaxis_title="Time",
        yaxis_title="Patterns Reviewed",
        height=300,
        margin=dict(l=20, r=20, t=40, b=20)
    )
    
    return fig_pie, review_velocity

# ============================================================================
# INITIALIZE SESSION STATE
# ============================================================================

# Core state
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

# API and cache
if 'api_key' not in st.session_state:
    st.session_state.api_key = ""
if 'viz_cache' not in st.session_state:
    st.session_state.viz_cache = {}

# Filter and results
if 'filter_results' not in st.session_state:
    st.session_state.filter_results = None
if 'uploaded_file_name' not in st.session_state:
    st.session_state.uploaded_file_name = None

# Search and filter
if 'search_query' not in st.session_state:
    st.session_state.search_query = ""
if 'active_filters' not in st.session_state:
    st.session_state.active_filters = []

# Session tracking
if 'session_start_time' not in st.session_state:
    st.session_state.session_start_time = datetime.now()
if 'session_name' not in st.session_state:
    st.session_state.session_name = f"Session_{datetime.now().strftime('%Y%m%d_%H%M%S')}"

# Bulk selection
if 'bulk_selected' not in st.session_state:
    st.session_state.bulk_selected = []

# Quick test state
if 'quick_test_results' not in st.session_state:
    st.session_state.quick_test_results = None
if 'quick_page_idx' not in st.session_state:
    st.session_state.quick_page_idx = 0

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

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
    
    # Apply search
    if search_query:
        mask = (
            filtered_df['SMARTS'].str.contains(search_query, case=False, na=False) |
            (filtered_df['Description'].str.contains(search_query, case=False, na=False) if 'Description' in filtered_df.columns else False)
        )
        filtered_df = filtered_df[mask]
    
    # Apply filters
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

# Navigation
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
    # Session duration
    duration = datetime.now() - st.session_state.session_start_time
    st.markdown(f"**Session Time:** {str(duration).split('.')[0]}")

st.divider()

# Settings expander
with st.expander("⚙️ Settings & Configuration", expanded=False):
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
            if st.button("💾 Save Session", use_container_width=True):
                session_json, filename = save_session()
                st.download_button(
                    "📥 Download Session",
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
                    st.success("Session loaded!")
                    st.rerun()
    
    with col_s3:
        st.subheader("Debug & Info")
        st.session_state.debug_mode = st.checkbox(
            "Debug mode",
            value=False,
            help="Show detailed error messages"
        )
        
        st.metric("Patterns Reviewed", len(st.session_state.decisions))
        st.metric("Cache Size", len(st.session_state.viz_cache))

st.divider()

# ============================================================================
# MODE: DASHBOARD
# ============================================================================

if mode == "Dashboard":
    st.header("📊 Dashboard Overview")
    
    # Quick stats
    col1, col2, col3, col4 = st.columns(4)
    
    stats = calculate_session_stats()
    
    with col1:
        total_reviewed = stats['total_reviewed'] if stats else 0
        st.markdown(metric_card("Total Reviewed", total_reviewed), unsafe_allow_html=True)
    
    with col2:
        pass_rate = stats['pass_rate'] if stats else 0
        st.markdown(metric_card("Pass Rate", f"{pass_rate:.1f}%"), unsafe_allow_html=True)
    
    with col3:
        block_rate = stats['block_rate'] if stats else 0
        st.markdown(metric_card("Block Rate", f"{block_rate:.1f}%"), unsafe_allow_html=True)
    
    with col4:
        total_patterns = len(st.session_state.smarts_data) if st.session_state.smarts_data is not None else 0
        st.markdown(metric_card("Total Patterns", total_patterns), unsafe_allow_html=True)
    
    st.write("")
    
    # Charts
    if stats:
        col_chart1, col_chart2 = st.columns(2)
        
        with col_chart1:
            fig_pie, _ = generate_analytics_charts()
            if fig_pie:
                st.plotly_chart(fig_pie, use_container_width=True)
        
        with col_chart2:
            _, fig_velocity = generate_analytics_charts()
            if fig_velocity:
                st.plotly_chart(fig_velocity, use_container_width=True)
    
    st.write("")
    
    # Quick actions
    st.subheader("Quick Actions")
    
    col_action1, col_action2, col_action3 = st.columns(3)
    
    with col_action1:
        st.markdown("""
        <div class="dashboard-card">
            <div class="dashboard-card-title">🔬 Pattern Curator</div>
            <div class="dashboard-card-subtitle">Review and classify SMARTS patterns</div>
        </div>
        """, unsafe_allow_html=True)
        if st.button("Start Curation", use_container_width=True, type="primary"):
            st.session_state.mode = "Pattern Curator"
            st.rerun()
    
    with col_action2:
        st.markdown("""
        <div class="dashboard-card">
            <div class="dashboard-card-title">⚡ Quick Tester</div>
            <div class="dashboard-card-subtitle">Test SMARTS patterns instantly</div>
        </div>
        """, unsafe_allow_html=True)
        if st.button("Quick Test", use_container_width=True):
            st.session_state.mode = "Quick Tester"
            st.rerun()
    
    with col_action3:
        st.markdown("""
        <div class="dashboard-card">
            <div class="dashboard-card-title">📈 Analytics</div>
            <div class="dashboard-card-subtitle">View detailed insights</div>
        </div>
        """, unsafe_allow_html=True)
        if st.button("View Analytics", use_container_width=True):
            st.session_state.mode = "Analytics"
            st.rerun()
    
    st.write("")
    
    # Recent activity
    if st.session_state.decisions and st.session_state.smarts_data is not None:
        st.subheader("Recent Activity")
        
        recent_decisions = list(st.session_state.decisions.items())[-5:]
        
        for idx, decision in reversed(recent_decisions):
            col_recent1, col_recent2 = st.columns([3, 1])
            
            with col_recent1:
                pattern = st.session_state.smarts_data.iloc[idx]['SMARTS']
                st.code(pattern[:80] + "..." if len(pattern) > 80 else pattern, language='text')
            
            with col_recent2:
                st.markdown(status_badge(decision), unsafe_allow_html=True)

# ============================================================================
# MODE: PATTERN CURATOR (Enhanced Visualizer)
# ============================================================================

elif mode == "Pattern Curator":
    st.header("🔬 Pattern Curator")
    
    # Search and filter bar
    col_search, col_filter, col_bulk = st.columns([3, 2, 1])
    
    with col_search:
        search_query = st.text_input(
            "🔍 Search patterns",
            placeholder="Search by SMARTS or description...",
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
    
    with col_bulk:
        if st.button("🗂️ Bulk Actions", use_container_width=True):
            st.info("Select patterns below to enable bulk operations")
    
    st.write("")
    
    # File upload
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
                st.success(f"✅ Loaded {len(df)} patterns")
            except Exception as e:
                st.error(f"Parse error: {str(e)}")
                st.stop()
        
        df = st.session_state.smarts_data
        
        # Apply filters
        filtered_df = filter_patterns(df, search_query, active_filters)
        
        if filtered_df is None or len(filtered_df) == 0:
            st.warning("No patterns match your search/filter criteria")
            st.stop()
        
        # Get current pattern from filtered set
        filtered_indices = filtered_df.index.tolist()
        
        # Find current position in filtered set
        if st.session_state.current_idx in filtered_indices:
            current_filtered_pos = filtered_indices.index(st.session_state.current_idx)
        else:
            current_filtered_pos = 0
            st.session_state.current_idx = filtered_indices[0]
        
        total = len(filtered_df)
        current = current_filtered_pos
        actual_idx = filtered_indices[current]
        
        # Progress bar
        progress = len(st.session_state.decisions) / len(df) if len(df) > 0 else 0
        st.markdown(progress_bar(len(st.session_state.decisions), len(df)), unsafe_allow_html=True)
        
        st.write("")
        
        # Stats row
        col1, col2, col3, col4, col5 = st.columns(5)
        with col1:
            st.metric("Reviewed", f"{len(st.session_state.decisions)}/{len(df)}")
        with col2:
            pass_count = sum(1 for d in st.session_state.decisions.values() if d == "OK")
            st.metric("Pass", pass_count)
            st.markdown(status_badge("OK"), unsafe_allow_html=True)
        with col3:
            yellow_count = sum(1 for d in st.session_state.decisions.values() if d == "Yellow")
            st.metric("Review", yellow_count)
            st.markdown(status_badge("Yellow"), unsafe_allow_html=True)
        with col4:
            amber_count = sum(1 for d in st.session_state.decisions.values() if d == "Amber")
            st.metric("Concern", amber_count)
            st.markdown(status_badge("Amber"), unsafe_allow_html=True)
        with col5:
            red_count = sum(1 for d in st.session_state.decisions.values() if d == "Red")
            st.metric("Block", red_count)
            st.markdown(status_badge("Red"), unsafe_allow_html=True)
        
        st.write("")
        
        if current < total:
            col_main, col_side = st.columns([2.5, 1])
            
            with col_main:
                smarts_pattern = df.iloc[actual_idx]['SMARTS']
                
                # Pattern info card
                st.markdown(f"""
                <div class="dashboard-card">
                    <div class="dashboard-card-title">Pattern {current + 1} of {total}</div>
                    <div class="dashboard-card-subtitle">Index: {actual_idx}</div>
                </div>
                """, unsafe_allow_html=True)
                
                st.code(smarts_pattern, language='text')
                
                if 'Description' in df.columns and pd.notna(df.iloc[actual_idx]['Description']):
                    st.caption(df.iloc[actual_idx]['Description'])
                
                # Visualization
                try:
                    pattern = Chem.MolFromSmarts(smarts_pattern)
                    if pattern is None:
                        st.error("⚠️ Invalid SMARTS pattern")
                    else:
                        image_displayed = False
                        if use_smartsplus and st.session_state.api_key:
                            with st.spinner("Loading SMARTS.plus visualization..."):
                                smartsplus_svg = get_smartsplus_image(smarts_pattern, st.session_state.api_key)
                                if smartsplus_svg:
                                    import re
                                    svg_fixed = re.sub(r'\s*width="[^"]*"', '', smartsplus_svg)
                                    svg_fixed = re.sub(r'\s*height="[^"]*"', '', svg_fixed)
                                    
                                    scaled_svg = f"""
                                    <div style="width: 800px; max-width: 100%;">
                                        {svg_fixed}
                                    </div>
                                    """
                                    st.markdown(scaled_svg, unsafe_allow_html=True)
                                    st.caption("🎨 SMARTS.plus")
                                    image_displayed = True
                        
                        if not image_displayed:
                            img = Draw.MolToImage(pattern, size=(400, 320))
                            st.image(img, width=400)
                except Exception as e:
                    st.error(f"Error: {str(e)}")
            
            with col_side:
                # Current decision
                if actual_idx in st.session_state.decisions:
                    current_decision = st.session_state.decisions[actual_idx]
                    st.markdown(f"**Current:** {status_badge(current_decision)}", unsafe_allow_html=True)
                    st.write("")
                
                # Decision buttons
                st.markdown("**Classify Pattern:**")
                
                if st.button("✓ PASS", use_container_width=True, key="pass", type="primary"):
                    st.session_state.decisions[actual_idx] = "OK"
                    # Auto-advance
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
                
                # Navigation
                col_n1, col_n2 = st.columns(2)
                with col_n1:
                    if st.button("⬅️ Prev", disabled=(current == 0), use_container_width=True):
                        st.session_state.current_idx = filtered_indices[current - 1]
                        st.rerun()
                with col_n2:
                    if st.button("Next ➡️", disabled=(current >= total - 1), use_container_width=True):
                        st.session_state.current_idx = filtered_indices[current + 1]
                        st.rerun()
                
                # Jump to
                jump_to = st.number_input(
                    "Jump to #", 
                    1, 
                    total, 
                    current + 1, 
                    key='jump',
                    label_visibility="collapsed"
                )
                if st.button("🎯 Jump", use_container_width=True):
                    st.session_state.current_idx = filtered_indices[jump_to - 1]
                    st.rerun()
        
        else:
            st.success("🎉 Review complete for filtered patterns!")
        
        # Export section
        if len(st.session_state.decisions) > 0:
            st.write("")
            st.divider()
            st.subheader("📥 Export Results")
            
            results_df = df.copy()
            results_df['Decision'] = results_df.index.map(
                lambda x: st.session_state.decisions.get(x, "NOT_REVIEWED")
            )
            
            col_e1, col_e2, col_e3, col_e4 = st.columns(4)
            
            with col_e1:
                st.download_button(
                    "📥 All Results",
                    results_df.to_csv(index=False),
                    "all_results.csv",
                    "text/csv",
                    use_container_width=True
                )
            
            with col_e2:
                concerning = results_df[results_df['Decision'].isin(['Amber', 'Red'])]
                if len(concerning) > 0:
                    st.download_button(
                        "⚠️ Concerning",
                        concerning.to_csv(index=False),
                        "concerning.csv",
                        "text/csv",
                        use_container_width=True
                    )
            
            with col_e3:
                ok_df = results_df[results_df['Decision'] == 'OK']
                if len(ok_df) > 0:
                    st.download_button(
                        "✅ Passed",
                        ok_df.to_csv(index=False),
                        "passed.csv",
                        "text/csv",
                        use_container_width=True
                    )
            
            with col_e4:
                # Export session
                session_json, filename = save_session()
                st.download_button(
                    "💾 Session",
                    session_json,
                    filename,
                    "application/json",
                    use_container_width=True
                )
    
    else:
        st.info("👆 Upload a CSV file with SMARTS patterns to begin curation")

# ============================================================================
# MODE: BATCH VALIDATOR
# ============================================================================

elif mode == "Batch Validator":
    st.header("🔬 Batch Validator")
    
    # Sub-mode selector
    validator_mode = st.radio(
        "Validation Type:",
        ["Batch Validation", "Quick Test"],
        horizontal=True,
        key="validator_submode"
    )
    
    st.write("")
    
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
                st.success(f"✅ Loaded {len(df)} SMARTS patterns")
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
            
            # Progress
            progress = current / total if total > 0 else 0
            st.markdown(progress_bar(current, total), unsafe_allow_html=True)
            
            st.write("")
            
            # Stats
            col_m1, col_m2 = st.columns(2)
            with col_m1:
                flagged = sum(1 for d in st.session_state.decisions.values() if d == "FLAGGED")
                st.metric("Flagged Patterns", flagged)
            with col_m2:
                checked = sum(1 for d in st.session_state.decisions.values() if d == "CHECKED")
                st.metric("Checked Patterns", checked)
            
            st.write("")
            
            if current < total:
                smarts_pattern = df.iloc[current]['SMARTS']
                
                col_main, col_side = st.columns([2.5, 1])
                
                with col_main:
                    st.markdown(f"""
                    <div class="dashboard-card">
                        <div class="dashboard-card-title">Pattern {current + 1} of {total}</div>
                    </div>
                    """, unsafe_allow_html=True)
                    
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
                            
                            # Match stats
                            match_rate = len(matches)/len(mol_df)*100
                            st.markdown(f"**{len(matches)}/{len(mol_df)} matches ({match_rate:.1f}%)**")
                            
                            # Match rate indicator
                            if match_rate > 50:
                                st.warning(f"⚠️ High match rate - pattern may be too broad")
                            elif match_rate == 0:
                                st.info("ℹ️ No matches - pattern may be too specific")
                            
                            if len(matches) > 0:
                                st.write("**Sample Matches:**")
                                cols = st.columns(3)
                                for i, (idx, mol) in enumerate(matches[:6]):
                                    with cols[i % 3]:
                                        img = Draw.MolToImage(mol, size=(150, 150), 
                                                             highlightAtoms=mol.GetSubstructMatch(pattern))
                                        st.image(img, width=150)
                                        if 'Name' in mol_df.columns:
                                            st.caption(mol_df.iloc[idx]['Name'], unsafe_allow_html=True)
                            else:
                                st.info("No matches found")
                    except Exception as e:
                        st.error(f"Error: {str(e)}")
                
                with col_side:
                    if current in st.session_state.decisions:
                        status = st.session_state.decisions[current]
                        if status == "FLAGGED":
                            st.markdown(status_badge("FLAGGED"), unsafe_allow_html=True)
                        else:
                            st.markdown(status_badge("CHECKED"), unsafe_allow_html=True)
                    
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
                        if st.button("⬅️", disabled=(current == 0), use_container_width=True, key='prev_v'):
                            st.session_state.current_idx -= 1
                            st.rerun()
                    with col_n2:
                        if st.button("➡️", disabled=(current >= total - 1), use_container_width=True, key='next_v'):
                            st.session_state.current_idx += 1
                            st.rerun()
            
            else:
                st.success("🎉 Validation complete!")
            
            if len(st.session_state.decisions) > 0:
                st.write("")
                st.divider()
                st.subheader("📥 Export Results")
                
                results_df = df.copy()
                results_df['Status'] = results_df.index.map(lambda x: st.session_state.decisions.get(x, "NOT_REVIEWED"))
                flagged = results_df[results_df['Status'] == 'FLAGGED']
                
                col_e1, col_e2 = st.columns(2)
                with col_e1:
                    st.download_button("📥 All Results", results_df.to_csv(index=False), "validation.csv", "text/csv", use_container_width=True)
                with col_e2:
                    if len(flagged) > 0:
                        st.download_button("🚩 Flagged Only", flagged.to_csv(index=False), "flagged.csv", "text/csv", use_container_width=True)
        
        else:
            st.info("👆 Upload both SMARTS patterns and molecules to begin validation")

# ============================================================================
# MODE: QUICK TESTER (Enhanced)
# ============================================================================

elif mode == "Quick Tester":
    st.header("⚡ Quick Tester")
    
    col_input, col_output = st.columns([1, 1])
    
    with col_input:
        st.subheader("Input")
        
        # Tabbed interface
        smarts_tab, library_tab = st.tabs(["Custom SMARTS", "Pattern Library"])
        
        with smarts_tab:
            quick_smarts = st.text_area(
                "SMARTS Pattern:",
                height=100,
                placeholder="[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1",
                help="Paste your SMARTS pattern to test",
                key="quick_smarts_custom"
            )
        
        with library_tab:
            st.write("**Common Patterns:**")
            
            pattern_library = {
                "Benzene Ring": "c1ccccc1",
                "Primary Amine": "[NX3;H2;!$(NC=O)]",
                "Carboxylic Acid": "C(=O)[OH]",
                "Aromatic Ring": "a",
                "Aliphatic Chain": "CCCC",
                "Ester": "C(=O)O",
                "Amide": "C(=O)N",
                "Ketone": "[CX3](=O)[#6]",
                "Aldehyde": "[CX3H1](=O)[#6]",
                "Ether": "[OD2]([#6])[#6]"
            }
            
            selected_pattern = st.selectbox(
                "Select a pattern",
                list(pattern_library.keys()),
                key="pattern_library_select"
            )
            
            if st.button("Use This Pattern", use_container_width=True):
                quick_smarts = pattern_library[selected_pattern]
                st.session_state.quick_smarts_custom = quick_smarts
                st.rerun()
        
        # Validate SMARTS in real-time
        if quick_smarts:
            pattern = Chem.MolFromSmarts(quick_smarts)
            if pattern:
                st.success("✓ Valid SMARTS pattern")
            else:
                st.error("✗ Invalid SMARTS pattern")
        
        st.divider()
        
        # Molecule source
        mol_source = st.radio(
            "Test Against:",
            ["Single SMILES", "Molecule File"],
            horizontal=True,
            key="quick_mol_source"
        )
        
        if mol_source == "Single SMILES":
            quick_smiles = st.text_input(
                "SMILES:",
                placeholder="CCO",
                help="Paste a single SMILES string",
                key="quick_smiles"
            )
            
            if quick_smarts and quick_smiles:
                if st.button("🚀 Test Match", type="primary", use_container_width=True):
                    try:
                        pattern = Chem.MolFromSmarts(quick_smarts)
                        mol = Chem.MolFromSmiles(quick_smiles)
                        
                        if pattern is None:
                            st.error("❌ Invalid SMARTS pattern")
                        elif mol is None:
                            st.error("❌ Invalid SMILES")
                        else:
                            is_match = mol.HasSubstructMatch(pattern)
                            
                            st.session_state.quick_test_results = {
                                'is_match': is_match,
                                'mol': mol,
                                'pattern': pattern,
                                'smarts': quick_smarts,
                                'smiles': quick_smiles
                            }
                    except Exception as e:
                        st.error(f"Error: {str(e)}")
        
        else:  # Molecule File
            quick_mol_file = st.file_uploader(
                "📁 Molecules (CSV/SDF)",
                type=['csv', 'sdf'],
                key='quick_mols',
                label_visibility="collapsed"
            )
            
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
            
            if st.button("🚀 Run Test", type="primary", use_container_width=True, disabled=not (quick_smarts and quick_mol_df is not None)):
                if quick_smarts and quick_mol_df is not None:
                    with st.spinner("Testing SMARTS pattern..."):
                        try:
                            pattern = Chem.MolFromSmarts(quick_smarts)
                            if pattern is None:
                                st.error("❌ Invalid SMARTS pattern")
                            else:
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
                        except Exception as e:
                            st.error(f"Error: {str(e)}")
    
    with col_output:
        st.subheader("Results")
        
        if st.session_state.quick_test_results:
            results = st.session_state.quick_test_results
            
            # Single SMILES result
            if 'is_match' in results:
                is_match = results['is_match']
                mol = results['mol']
                pattern = results['pattern']
                
                if is_match:
                    st.markdown("""
                    <div class="toast-success">
                        <strong>✅ MATCH FOUND</strong><br>
                        The SMARTS pattern matches this molecule
                    </div>
                    """, unsafe_allow_html=True)
                    
                    match_atoms = mol.GetSubstructMatch(pattern)
                    img = Draw.MolToImage(mol, size=(400, 400), highlightAtoms=match_atoms)
                else:
                    st.markdown("""
                    <div class="toast-warning">
                        <strong>❌ NO MATCH</strong><br>
                        The SMARTS pattern does not match this molecule
                    </div>
                    """, unsafe_allow_html=True)
                    
                    img = Draw.MolToImage(mol, size=(400, 400))
                
                st.image(img, use_column_width=True)
            
            # File result
            elif 'matches' in results:
                matches = results['matches']
                total = results['total']
                
                # Summary metrics
                col_r1, col_r2, col_r3 = st.columns(3)
                with col_r1:
                    st.metric("Total Molecules", total)
                with col_r2:
                    st.metric("Matches", len(matches))
                with col_r3:
                    match_rate = len(matches)/total*100 if total > 0 else 0
                    st.metric("Match Rate", f"{match_rate:.1f}%")
                
                st.write("")
                
                if len(matches) > 0:
                    # Pagination
                    per_page = 6
                    total_pages = (len(matches) + per_page - 1) // per_page
                    current_page = st.session_state.quick_page_idx
                    
                    start_idx = current_page * per_page
                    end_idx = min(start_idx + per_page, len(matches))
                    page_matches = matches[start_idx:end_idx]
                    
                    st.write(f"**Showing matches {start_idx + 1}-{end_idx} of {len(matches)}**")
                    
                    cols = st.columns(3)
                    for i, (idx, mol, row) in enumerate(page_matches):
                        with cols[i % 3]:
                            match_atoms = mol.GetSubstructMatch(results['pattern'])
                            img = Draw.MolToImage(mol, size=(150, 150), highlightAtoms=match_atoms)
                            st.image(img, width=150)
                            
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
                            st.write(f"**Page {current_page + 1} of {total_pages}**")
                        
                        with col_p3:
                            if st.button("Next ➡️", disabled=(current_page >= total_pages - 1), use_container_width=True):
                                st.session_state.quick_page_idx += 1
                                st.rerun()
                    
                    # Export
                    st.write("")
                    match_indices = [idx for idx, _, _ in matches]
                    export_df = results['mol_df'].iloc[match_indices].copy()
                    
                    export_cols = [col for col in export_df.columns if col not in ['Mol', 'ROMol']]
                    
                    if 'SMILES' not in export_cols:
                        export_df['SMILES'] = export_df['Mol'].apply(lambda m: Chem.MolToSmiles(m) if m else '')
                        export_cols = ['SMILES'] + export_cols
                    
                    st.download_button(
                        "📥 Download Matches",
                        export_df[export_cols].to_csv(index=False),
                        "smarts_matches.csv",
                        "text/csv",
                        use_container_width=True,
                        type="primary"
                    )
                else:
                    st.info("No matches found for this SMARTS pattern")
        else:
            st.info("👈 Configure your test and click the Run button to see results")

# ============================================================================
# MODE: PRODUCTION FILTER (Gen AI Filter)
# ============================================================================

elif mode == "Production Filter":
    st.header("🚀 Production Filter")
    st.caption("Apply approved SMARTS patterns to filter Gen AI output")
    
    col_u1, col_u2 = st.columns(2)
    with col_u1:
        filter_smarts_file = st.file_uploader(
            "📁 Approved SMARTS (CSV from Pattern Curator)",
            type=['csv'],
            key='filter_smarts',
            help="Upload CSV with 'Decision' column"
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
        "Apply patterns with these classifications:",
        ["Red", "Amber", "Yellow", "OK"],
        default=["Red", "Amber"],
        help="Which classifications to use as filters"
    )
    
    # Load filter patterns
    filter_patterns = None
    if filter_smarts_file:
        try:
            filter_df = pd.read_csv(filter_smarts_file)
            
            if 'Decision' not in filter_df.columns:
                st.warning("No 'Decision' column found. Using all patterns.")
                filter_patterns = filter_df
            else:
                filter_patterns = filter_df[filter_df['Decision'].isin(filter_severity)]
                st.success(f"✅ Loaded {len(filter_patterns)} filter patterns ({', '.join(filter_severity)})")
                
        except Exception as e:
            st.error(f"Error loading SMARTS: {str(e)}")
    
    # Load molecules
    gen_ai_mols = None
    if gen_ai_file:
        try:
            if gen_ai_file.name.endswith('.sdf'):
                gen_ai_mols = PandasTools.LoadSDF(gen_ai_file, molColName='Mol')
                st.success(f"✅ Loaded {len(gen_ai_mols)} molecules from SDF")
            else:
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
        
        if st.button("🚀 Run Production Filter", type="primary", use_container_width=True):
            
            with st.spinner("Filtering molecules..."):
                
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
                
                passed = []
                failed = []
                rejection_reasons = []
                
                progress_bar_placeholder = st.empty()
                
                for idx, row in gen_ai_mols.iterrows():
                    mol = row.get('Mol') or row.get('ROMol')
                    if mol is None:
                        continue
                    
                    # Update progress
                    progress = (idx + 1) / len(gen_ai_mols)
                    progress_bar_placeholder.progress(progress, text=f"Processing molecule {idx + 1}/{len(gen_ai_mols)}")
                    
                    caught_by = []
                    for p in compiled_patterns:
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
                            'descriptions': [p['description'] for p in caught_by],
                            'decisions': [p['decision'] for p in caught_by]
                        })
                
                progress_bar_placeholder.empty()
                
                st.session_state.filter_results = {
                    'passed': passed,
                    'failed': failed,
                    'rejection_reasons': rejection_reasons,
                    'total': len(gen_ai_mols),
                    'patterns_used': len(compiled_patterns)
                }
                
                st.rerun()
    
    # Display results
    if st.session_state.filter_results:
        results = st.session_state.filter_results
        
        st.write("")
        st.divider()
        st.subheader("📊 Filtering Results")
        
        # Summary metrics with gradient cards
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.markdown(metric_card("Total Molecules", results['total']), unsafe_allow_html=True)
        
        with col2:
            pass_rate = len(results['passed'])/results['total']*100
            st.markdown(metric_card("Passed", f"{len(results['passed'])} ({pass_rate:.1f}%)"), unsafe_allow_html=True)
        
        with col3:
            fail_rate = len(results['failed'])/results['total']*100
            st.markdown(metric_card("Failed", f"{len(results['failed'])} ({fail_rate:.1f}%)"), unsafe_allow_html=True)
        
        with col4:
            st.markdown(metric_card("Filters Applied", results['patterns_used']), unsafe_allow_html=True)
        
        st.write("")
        
        # Export buttons
        col_e1, col_e2 = st.columns(2)
        
        with col_e1:
            if len(results['passed']) > 0:
                passed_mols = gen_ai_mols.iloc[results['passed']].copy()
                
                if 'SMILES' in passed_mols.columns:
                    export_df = passed_mols[['SMILES']].copy()
                else:
                    export_df = pd.DataFrame({
                        'SMILES': passed_mols['Mol'].apply(lambda m: Chem.MolToSmiles(m) if m else '')
                    })
                
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
            if len(results['failed']) > 0:
                rejection_df = pd.DataFrame(results['rejection_reasons'])
                
                st.download_button(
                    "📥 Download Rejection Report",
                    rejection_df.to_csv(index=False),
                    "rejection_report.csv",
                    "text/csv",
                    use_container_width=True
                )
        
        # Sample rejections
        if len(results['failed']) > 0:
            with st.expander(f"🔍 View Sample Rejections (first 5 of {len(results['failed'])})"):
                for i, reason in enumerate(results['rejection_reasons'][:5]):
                    st.markdown(f"""
                    <div class="dashboard-card">
                        <div class="dashboard-card-title">Molecule {reason['molecule_idx']}</div>
                        <div class="dashboard-card-subtitle">Caught by {reason['num_violations']} pattern(s)</div>
                    </div>
                    """, unsafe_allow_html=True)
                    
                    for j, (pattern, desc, decision) in enumerate(zip(reason['patterns'], reason['descriptions'], reason['decisions'])):
                        col_rej1, col_rej2 = st.columns([3, 1])
                        with col_rej1:
                            st.code(pattern, language='text')
                            if desc:
                                st.caption(desc)
                        with col_rej2:
                            st.markdown(status_badge(decision), unsafe_allow_html=True)
                    
                    if i < 4:
                        st.divider()
    
    elif filter_patterns is not None and gen_ai_mols is not None and len(filter_patterns) == 0:
        st.warning("⚠️ No patterns selected. Choose at least one classification category to filter.")
    
    else:
        st.info("👆 Upload approved SMARTS patterns and Gen AI molecules, then click 'Run Production Filter'")

# ============================================================================
# MODE: ANALYTICS
# ============================================================================

elif mode == "Analytics":
    st.header("📈 Analytics Dashboard")
    
    stats = calculate_session_stats()
    
    if stats:
        # Key metrics
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.markdown(metric_card("Total Reviewed", stats['total_reviewed']), unsafe_allow_html=True)
        
        with col2:
            st.markdown(metric_card("Pass Rate", f"{stats['pass_rate']:.1f}%"), unsafe_allow_html=True)
        
        with col3:
            st.markdown(metric_card("Block Rate", f"{stats['block_rate']:.1f}%"), unsafe_allow_html=True)
        
        with col4:
            session_duration = datetime.now() - st.session_state.session_start_time
            minutes = int(session_duration.total_seconds() / 60)
            st.markdown(metric_card("Session Duration", f"{minutes} min"), unsafe_allow_html=True)
        
        st.write("")
        
        # Charts
        col_chart1, col_chart2 = st.columns(2)
        
        with col_chart1:
            fig_pie, _ = generate_analytics_charts()
            if fig_pie:
                st.plotly_chart(fig_pie, use_container_width=True)
        
        with col_chart2:
            _, fig_velocity = generate_analytics_charts()
            if fig_velocity:
                st.plotly_chart(fig_velocity, use_container_width=True)
        
        st.write("")
        st.divider()
        
        # Detailed breakdown
        st.subheader("Decision Breakdown")
        
        col_detail1, col_detail2, col_detail3, col_detail4 = st.columns(4)
        
        decisions_count = stats['decisions_count']
        
        with col_detail1:
            pass_count = decisions_count.get('OK', 0)
            st.metric("Pass", pass_count)
            st.markdown(status_badge("OK"), unsafe_allow_html=True)
            if stats['total_reviewed'] > 0:
                st.caption(f"{pass_count/stats['total_reviewed']*100:.1f}% of total")
        
        with col_detail2:
            yellow_count = decisions_count.get('Yellow', 0)
            st.metric("Review", yellow_count)
            st.markdown(status_badge("Yellow"), unsafe_allow_html=True)
            if stats['total_reviewed'] > 0:
                st.caption(f"{yellow_count/stats['total_reviewed']*100:.1f}% of total")
        
        with col_detail3:
            amber_count = decisions_count.get('Amber', 0)
            st.metric("Concern", amber_count)
            st.markdown(status_badge("Amber"), unsafe_allow_html=True)
            if stats['total_reviewed'] > 0:
                st.caption(f"{amber_count/stats['total_reviewed']*100:.1f}% of total")
        
        with col_detail4:
            red_count = decisions_count.get('Red', 0)
            st.metric("Block", red_count)
            st.markdown(status_badge("Red"), unsafe_allow_html=True)
            if stats['total_reviewed'] > 0:
                st.caption(f"{red_count/stats['total_reviewed']*100:.1f}% of total")
        
        st.write("")
        st.divider()
        
        # Recommendations
        st.subheader("🎯 Insights & Recommendations")
        
        if stats['block_rate'] > 20:
            st.markdown("""
            <div class="toast-warning">
                <strong>⚠️ High Block Rate</strong><br>
                Over 20% of patterns are being blocked. Consider reviewing your filter criteria to ensure they're not too restrictive.
            </div>
            """, unsafe_allow_html=True)
        
        if stats['pass_rate'] > 80:
            st.markdown("""
            <div class="toast-success">
                <strong>✅ High Pass Rate</strong><br>
                Great! Over 80% of patterns are passing. Your curation process is working well.
            </div>
            """, unsafe_allow_html=True)
        
        if stats['total_reviewed'] < 50:
            st.info("💡 **Tip:** Review at least 50 patterns for meaningful statistics.")
    
    else:
        st.info("No analytics data available yet. Start reviewing patterns in Pattern Curator mode!")

# ============================================================================
# FOOTER
# ============================================================================

st.write("")
st.divider()
st.caption("🧬 SMARTS Toolkit Pro | Built for Digital Chemistry Teams")
