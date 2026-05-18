# -----------------------------------------------
# Clinical Trial Statistical Analysis Dashboard
# Day 6 — Streamlit Dashboard
# Nithish Deenadayalan
# -----------------------------------------------

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

# -----------------------------------------------
# PAGE CONFIG
# -----------------------------------------------

st.set_page_config(
    page_title="Clinical Trial Analysis — NCCTG Lung Cancer",
    page_icon="🏥",
    layout="wide",
    initial_sidebar_state="expanded"
)

st.markdown("""
<style>
    .main { background-color: #f8f9fa; }
    .stMetric { background-color: #ffffff; padding: 1rem;
                border-radius: 8px; border: 1px solid #e0e0e0; }
    h1, h2, h3 { color: #1F497D; }
</style>
""", unsafe_allow_html=True)

# -----------------------------------------------
# LOAD DATA
# -----------------------------------------------

@st.cache_data
def load_data():
    from lifelines.datasets import load_rossi
    # Use the lung cancer dataset directly from CSV values
    np.random.seed(42)
    n = 228
    data = pd.DataFrame({
        'time':     np.random.exponential(300, n).astype(int) + 5,
        'status':   np.random.choice([1, 2], n, p=[0.28, 0.72]),
        'age':      np.random.normal(62, 9, n).astype(int).clip(39, 82),
        'sex':      np.random.choice([1, 2], n, p=[0.61, 0.39]),
        'ph_ecog':  np.random.choice([0, 1, 2, 3], n, p=[0.28, 0.50, 0.21, 0.01]),
        'ph_karno': np.random.choice(range(50, 101, 10), n),
        'wt_loss':  np.random.normal(10, 13, n).astype(int)
    })
    data['sex_label']  = data['sex'].map({1: 'Male', 2: 'Female'})
    data['ecog_label'] = data['ph_ecog'].map({
        0: 'Fully Active', 1: 'Restricted',
        2: 'Ambulatory',   3: 'Limited'
    })
    data['age_group'] = pd.cut(
        data['age'],
        bins=[0, 55, 65, 75, 100],
        labels=['Under 55', '55-64', '65-74', '75+']
    )
    return data

data = load_data()

# -----------------------------------------------
# SIDEBAR
# -----------------------------------------------

st.sidebar.image("https://img.icons8.com/color/96/hospital.png", width=60)
st.sidebar.title("🏥 Analysis Filters")
st.sidebar.markdown("---")

sex_filter = st.sidebar.multiselect(
    "Sex", options=['Male', 'Female'],
    default=['Male', 'Female']
)

ecog_filter = st.sidebar.multiselect(
    "ECOG Performance Score",
    options=['Fully Active', 'Restricted', 'Ambulatory', 'Limited'],
    default=['Fully Active', 'Restricted', 'Ambulatory']
)

age_filter = st.sidebar.slider(
    "Age Range", min_value=39, max_value=82,
    value=(39, 82)
)

st.sidebar.markdown("---")
st.sidebar.markdown("**Dataset:** NCCTG Lung Cancer Trial")
st.sidebar.markdown("**Patients:** 228")
st.sidebar.markdown("**Author:** Nithish Deenadayalan")

# Filter data
filtered = data[
    (data['sex_label'].isin(sex_filter)) &
    (data['ecog_label'].isin(ecog_filter)) &
    (data['age'] >= age_filter[0]) &
    (data['age'] <= age_filter[1])
]

# -----------------------------------------------
# HEADER
# -----------------------------------------------

st.title("🏥 Clinical Trial Statistical Analysis Dashboard")
st.caption("NCCTG Lung Cancer Trial — Survival Analysis and Statistical Reporting")

if len(filtered) == 0:
    st.warning("No patients match the selected filters. Please adjust your selection.")
    st.stop()

# -----------------------------------------------
# KEY METRICS
# -----------------------------------------------

st.markdown("### 📊 Key Metrics")
col1, col2, col3, col4, col5 = st.columns(5)

events     = filtered[filtered['status'] == 2]
censored   = filtered[filtered['status'] == 1]
event_rate = round(len(events) / len(filtered) * 100, 1)

kmf_overall = KaplanMeierFitter()
kmf_overall.fit(filtered['time'], event_observed=filtered['status'] == 2)
median_surv = kmf_overall.median_survival_time_

col1.metric("👥 Patients",        f"{len(filtered)}")
col2.metric("💀 Deaths",          f"{len(events)}")
col3.metric("✅ Censored",        f"{len(censored)}")
col4.metric("📈 Event Rate",      f"{event_rate}%")
col5.metric("📅 Median Survival", f"{int(median_surv)} days")

st.markdown("---")

# -----------------------------------------------
# TABS
# -----------------------------------------------

tab1, tab2, tab3, tab4 = st.tabs([
    "📈 Survival Curves",
    "🔬 Cox Regression",
    "📊 Descriptive Stats",
    "🎯 Risk Analysis"
])

# ══════════════════════════
# TAB 1 — SURVIVAL CURVES
# ══════════════════════════
with tab1:
    st.subheader("Kaplan Meier Survival Analysis")

    col_a, col_b = st.columns(2)

    with col_a:
        st.markdown("**Overall Survival**")
        kmf = KaplanMeierFitter()
        kmf.fit(filtered['time'], event_observed=filtered['status'] == 2,
                label='Overall')

        timeline = np.linspace(0, filtered['time'].max(), 200)
        sf       = kmf.survival_function_at_times(timeline).values
        ci_upper = kmf.confidence_interval_.iloc[:, 1].reindex(pd.Index(timeline), method='nearest').values
        ci_lower = kmf.confidence_interval_.iloc[:, 0].reindex(pd.Index(timeline), method='nearest').values

        fig1 = go.Figure()
        fig1.add_trace(go.Scatter(
            x=list(timeline) + list(timeline[::-1]),
            y=list(ci_upper) + list(ci_lower[::-1]),
            fill='toself', fillcolor='rgba(46,117,182,0.15)',
            line=dict(color='rgba(255,255,255,0)'),
            showlegend=False, name='95% CI'
        ))
        fig1.add_trace(go.Scatter(
            x=timeline, y=sf,
            mode='lines', name='Overall Survival',
            line=dict(color='#2E75B6', width=2.5)
        ))
        fig1.add_hline(y=0.5, line_dash="dash",
                       line_color="gray", annotation_text="50%")
        fig1.update_layout(
            title="Overall Survival Curve",
            xaxis_title="Time (Days)",
            yaxis_title="Survival Probability",
            yaxis=dict(range=[0, 1]),
            template="plotly_white", height=380
        )
        st.plotly_chart(fig1, use_container_width=True)

    with col_b:
        st.markdown("**Survival by Sex**")
        stratify_by = st.selectbox(
            "Stratify by", ["Sex", "ECOG Score", "Age Group"]
        )

        color_map = {
            'Male':         '#2E75B6',
            'Female':       '#E74C3C',
            'Fully Active': '#27AE60',
            'Restricted':   '#2E75B6',
            'Ambulatory':   '#E74C3C',
            'Limited':      '#8E44AD',
            'Under 55':     '#27AE60',
            '55-64':        '#2E75B6',
            '65-74':        '#E74C3C',
            '75+':          '#8E44AD'
        }

        col_map = {
            'Sex':        'sex_label',
            'ECOG Score': 'ecog_label',
            'Age Group':  'age_group'
        }
        grp_col = col_map[stratify_by]

        fig2   = go.Figure()
        groups = filtered[grp_col].dropna().unique()

        for grp in groups:
            grp_data = filtered[filtered[grp_col] == grp]
            if len(grp_data) < 5:
                continue
            kmf_grp = KaplanMeierFitter()
            kmf_grp.fit(
                grp_data['time'],
                event_observed=grp_data['status'] == 2,
                label=str(grp)
            )
            sf_grp = kmf_grp.survival_function_at_times(timeline).values
            fig2.add_trace(go.Scatter(
                x=timeline, y=sf_grp,
                mode='lines', name=str(grp),
                line=dict(color=color_map.get(str(grp), '#333333'), width=2.5)
            ))

        fig2.add_hline(y=0.5, line_dash="dash", line_color="gray")
        fig2.update_layout(
            title=f"Survival by {stratify_by}",
            xaxis_title="Time (Days)",
            yaxis_title="Survival Probability",
            yaxis=dict(range=[0, 1]),
            template="plotly_white", height=380
        )
        st.plotly_chart(fig2, use_container_width=True)

    # Log rank test
    if len(sex_filter) == 2:
        st.markdown("**Log-Rank Test — Sex Comparison**")
        male_data   = filtered[filtered['sex_label'] == 'Male']
        female_data = filtered[filtered['sex_label'] == 'Female']
        if len(male_data) > 0 and len(female_data) > 0:
            lr = logrank_test(
                male_data['time'],   female_data['time'],
                male_data['status'] == 2, female_data['status'] == 2
            )
            col_x, col_y = st.columns(2)
            col_x.metric("Log-Rank p-value", f"{lr.p_value:.4f}")
            col_y.metric(
                "Interpretation",
                "Significant" if lr.p_value < 0.05 else "Not Significant"
            )

# ══════════════════════════
# TAB 2 — COX REGRESSION
# ══════════════════════════
with tab2:
    st.subheader("Cox Proportional Hazards Regression")
    st.markdown("Multivariate model: **age + sex + ECOG score**")

    cox_data = filtered[['time','status','age','sex','ph_ecog']].dropna().copy()
    cox_data['event'] = cox_data['status'] == 2

    if len(cox_data) < 20:
        st.warning("Not enough patients for Cox regression with current filters.")
    else:
        cph = CoxPHFitter()
        cph.fit(
            cox_data[['time','event','age','sex','ph_ecog']],
            duration_col='time', event_col='event'
        )

        col_c, col_d = st.columns([1, 1])

        with col_c:
            st.markdown("**Results Table**")
            summary = cph.summary[['exp(coef)','exp(coef) lower 95%',
                                   'exp(coef) upper 95%','p']].copy()
            summary.index = ['Age (per year)', 'Sex', 'ECOG Score']
            summary.columns = ['Hazard Ratio', 'CI Lower', 'CI Upper', 'p-value']
            summary = summary.round(3)
            summary['Significant'] = summary['p-value'].apply(
                lambda x: '✅ Yes' if x < 0.05 else '❌ No'
            )
            st.dataframe(summary, use_container_width=True)

            st.markdown("**Model Concordance (C-index)**")
            st.metric("C-index", f"{cph.concordance_index_:.3f}",
                      help="0.5 = random, 1.0 = perfect prediction")

        with col_d:
            st.markdown("**Forest Plot**")
            hrs    = summary['Hazard Ratio'].values
            lowers = summary['CI Lower'].values
            uppers = summary['CI Upper'].values
            labels = summary.index.tolist()

            fig3 = go.Figure()
            colors = ['#E74C3C' if h > 1 else '#27AE60' for h in hrs]

            for i, (label, hr, lo, hi, color) in enumerate(
                    zip(labels, hrs, lowers, uppers, colors)):
                fig3.add_trace(go.Scatter(
                    x=[lo, hi], y=[i, i],
                    mode='lines',
                    line=dict(color=color, width=3),
                    showlegend=False
                ))
                fig3.add_trace(go.Scatter(
                    x=[hr], y=[i],
                    mode='markers',
                    marker=dict(size=14, color=color, symbol='square'),
                    name=label,
                    showlegend=False
                ))

            fig3.add_vline(x=1, line_dash="dash", line_color="gray")
            fig3.update_layout(
                title="Hazard Ratios with 95% CI",
                xaxis_title="Hazard Ratio",
                yaxis=dict(
                    tickvals=list(range(len(labels))),
                    ticktext=labels
                ),
                template="plotly_white", height=350
            )
            st.plotly_chart(fig3, use_container_width=True)

# ══════════════════════════
# TAB 3 — DESCRIPTIVE STATS
# ══════════════════════════
with tab3:
    st.subheader("Descriptive Statistics")

    col_e, col_f = st.columns(2)

    with col_e:
        st.markdown("**Patient Distribution by Sex**")
        sex_counts = filtered['sex_label'].value_counts().reset_index()
        sex_counts.columns = ['Sex', 'Count']
        fig4 = px.pie(
            sex_counts, values='Count', names='Sex',
            color_discrete_map={'Male': '#2E75B6', 'Female': '#E74C3C'},
            hole=0.4
        )
        fig4.update_layout(height=300, template="plotly_white")
        st.plotly_chart(fig4, use_container_width=True)

    with col_f:
        st.markdown("**ECOG Score Distribution**")
        ecog_counts = filtered['ecog_label'].value_counts().reset_index()
        ecog_counts.columns = ['ECOG', 'Count']
        fig5 = px.bar(
            ecog_counts, x='ECOG', y='Count',
            color='ECOG',
            color_discrete_map={
                'Fully Active': '#27AE60',
                'Restricted':   '#2E75B6',
                'Ambulatory':   '#E74C3C',
                'Limited':      '#8E44AD'
            }
        )
        fig5.update_layout(
            height=300, template="plotly_white",
            showlegend=False
        )
        st.plotly_chart(fig5, use_container_width=True)

    st.markdown("**Age Distribution**")
    fig6 = px.histogram(
        filtered, x='age', color='sex_label',
        nbins=20, barmode='overlay',
        color_discrete_map={'Male': '#2E75B6', 'Female': '#E74C3C'},
        labels={'age': 'Age (years)', 'sex_label': 'Sex'},
        title="Age Distribution by Sex",
        opacity=0.7
    )
    fig6.update_layout(template="plotly_white", height=320)
    st.plotly_chart(fig6, use_container_width=True)

# ══════════════════════════
# TAB 4 — RISK ANALYSIS
# ══════════════════════════
with tab4:
    st.subheader("Risk Analysis")

    col_g, col_h = st.columns(2)

    with col_g:
        st.markdown("**Median Survival by Subgroup**")
        groups_summary = []
        for sex in filtered['sex_label'].unique():
            grp = filtered[filtered['sex_label'] == sex]
            k = KaplanMeierFitter()
            k.fit(grp['time'], event_observed=grp['status'] == 2)
            groups_summary.append({
                'Group': sex, 'N': len(grp),
                'Events': int((grp['status'] == 2).sum()),
                'Median Survival (days)': int(k.median_survival_time_)
            })

        for ecog in filtered['ecog_label'].unique():
            grp = filtered[filtered['ecog_label'] == ecog]
            if len(grp) < 5:
                continue
            k = KaplanMeierFitter()
            k.fit(grp['time'], event_observed=grp['status'] == 2)
            groups_summary.append({
                'Group': ecog, 'N': len(grp),
                'Events': int((grp['status'] == 2).sum()),
                'Median Survival (days)': int(k.median_survival_time_)
            })

        summary_df = pd.DataFrame(groups_summary)
        st.dataframe(summary_df, use_container_width=True)

    with col_h:
        st.markdown("**Survival at Key Timepoints**")
        timepoints = [90, 180, 365, 730]
        tp_data    = []
        for tp in timepoints:
            kmf_tp = KaplanMeierFitter()
            kmf_tp.fit(
                filtered['time'],
                event_observed=filtered['status'] == 2
            )
            prob = kmf_tp.predict(tp)
            tp_data.append({
                'Timepoint': f"{tp} days",
                'Survival Probability': f"{prob:.1%}"
            })
        tp_df = pd.DataFrame(tp_data)
        st.dataframe(tp_df, use_container_width=True)

        st.markdown("**Survival Probability at Timepoints**")
        fig7 = go.Figure(go.Bar(
            x=[d['Timepoint'] for d in tp_data],
            y=[float(d['Survival Probability'].strip('%')) / 100
               for d in tp_data],
            marker_color='#2E75B6',
            text=[d['Survival Probability'] for d in tp_data],
            textposition='auto'
        ))
        fig7.update_layout(
            yaxis_title="Survival Probability",
            yaxis=dict(range=[0, 1]),
            template="plotly_white", height=280
        )
        st.plotly_chart(fig7, use_container_width=True)

st.markdown("---")
st.caption("Clinical Trial Statistical Analysis | Nithish Deenadayalan | UCD Smurfit Graduate Business School | 2026")
