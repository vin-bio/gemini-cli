
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge, Lasso
from sklearn.metrics import mean_absolute_error, r2_score
import xgboost as xgb
import lightgbm as lgb
import warnings

# Suppress future warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# Use a publication-ready style for plots
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_context("paper", font_scale=1.5)

# --- 1. Load and Preprocess Data ---
def load_and_preprocess_data(file_path):
    print("Loading and preprocessing data...")
    adata = anndata.read_h5ad(file_path)
    
    # Convert to DataFrame for easier manipulation
    df = adata.to_df()
    df['age'] = adata.obs['age'].values
    df['sex'] = adata.obs['sex'].values
    df['organ'] = adata.obs['organ'].values
    
    # Convert age to numeric and drop rows where age is not a valid number
    df['age'] = pd.to_numeric(df['age'], errors='coerce')
    df.dropna(subset=['age'], inplace=True)

    # For this regression task, we will focus on age prediction from gene expression
    # We will not use sex and organ as features for now
    X = df.drop(['age', 'sex', 'organ'], axis=1)
    y = df['age']
    
    # Scale the data
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    
    # PCA for dimensionality reduction
    pca = PCA(n_components=50)
    X_pca = pca.fit_transform(X_scaled)
    
    print("Data loaded and preprocessed.")
    return X_pca, y, pca, df

# --- 2. Train and Evaluate Models ---
def train_and_evaluate_models(X_train, X_test, y_train, y_test):
    models = {
        'Ridge': Ridge(),
        'Lasso': Lasso(),
        'XGBoost': xgb.XGBRegressor(objective='reg:squarederror', n_estimators=100, random_state=42),
        'LightGBM': lgb.LGBMRegressor(random_state=42)
    }
    
    results = {}
    print("Training and evaluating models...")
    for name, model in models.items():
        print(f"Training {name}...")
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        mae = mean_absolute_error(y_test, y_pred)
        r2 = r2_score(y_test, y_pred)
        results[name] = {'model': model, 'y_pred': y_pred, 'mae': mae, 'r2': r2}
    
    print("Model training and evaluation complete.")
    return results

# --- 3. Generate Refined Figures ---
def generate_refined_figures(results, pca, df, y_test):
    print("Generating refined figures...")
    
    # Define a color palette
    palette = sns.color_palette("viridis", 4)

    # Figure 1: Data Overview
    fig1, axes1 = plt.subplots(2, 2, figsize=(12, 10))
    fig1.suptitle('Figure 1: Dataset Overview and Preprocessing', fontsize=20, fontweight='bold')
    
    # Panel A: Age Distribution
    sns.histplot(df['age'], bins=30, ax=axes1[0, 0], kde=True, color=palette[0])
    axes1[0, 0].set_title('Age Distribution', fontsize=16)
    axes1[0, 0].set_xlabel('Age', fontsize=12)
    axes1[0, 0].set_ylabel('Frequency', fontsize=12)
    axes1[0, 0].text(-0.1, 1.1, 'A', transform=axes1[0, 0].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

    # Panel B: Sex Distribution
    sex_counts = df['sex'].value_counts()
    axes1[0, 1].pie(sex_counts, labels=sex_counts.index, autopct='%1.1f%%', startangle=90, colors=sns.color_palette("pastel"))
    axes1[0, 1].set_title('Sex Distribution', fontsize=16)
    axes1[0, 1].text(-0.1, 1.1, 'B', transform=axes1[0, 1].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

    # Panel C: Organ Distribution
    organ_counts = df['organ'].value_counts().nlargest(10)
    sns.barplot(x=organ_counts.values, y=organ_counts.index, ax=axes1[1, 0], palette="viridis")
    axes1[1, 0].set_title('Top 10 Organ Distribution', fontsize=16)
    axes1[1, 0].set_xlabel('Count', fontsize=12)
    axes1[1, 0].set_ylabel('Organ', fontsize=12)
    axes1[1, 0].text(-0.1, 1.1, 'C', transform=axes1[1, 0].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

    # Panel D: PCA Explained Variance
    axes1[1, 1].plot(np.cumsum(pca.explained_variance_ratio_), color=palette[3], lw=2)
    axes1[1, 1].set_xlabel('Number of Components', fontsize=12)
    axes1[1, 1].set_ylabel('Cumulative Explained Variance', fontsize=12)
    axes1[1, 1].set_title('PCA Explained Variance', fontsize=16)
    axes1[1, 1].grid(True, which='both', linestyle='--', linewidth=0.5)
    axes1[1, 1].text(-0.1, 1.1, 'D', transform=axes1[1, 1].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
    
    fig1.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig1.savefig('figure1_data_overview_refined.png', dpi=300)
    plt.close(fig1)

    # Figure 2: Classical Models
    fig2, axes2 = plt.subplots(2, 2, figsize=(12, 10))
    fig2.suptitle('Figure 2: Classical Regression Models Performance', fontsize=20, fontweight='bold')
    
    # Ridge
    axes2[0, 0].scatter(y_test, results['Ridge']['y_pred'], alpha=0.5, color=palette[0])
    axes2[0, 0].plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'r--', lw=2)
    axes2[0, 0].set_title(f"Ridge Regression (MAE: {results['Ridge']['mae']:.2f}, R2: {results['Ridge']['r2']:.2f})", fontsize=14)
    axes2[0, 0].set_xlabel('Actual Age', fontsize=12)
    axes2[0, 0].set_ylabel('Predicted Age', fontsize=12)
    axes2[0, 0].text(-0.1, 1.1, 'A', transform=axes2[0, 0].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

    sns.residplot(x=y_test, y=results['Ridge']['y_pred'], ax=axes2[0, 1], color=palette[0])
    axes2[0, 1].set_title('Ridge Residuals', fontsize=14)
    axes2[0, 1].set_xlabel('Actual Age', fontsize=12)
    axes2[0, 1].set_ylabel('Residuals', fontsize=12)
    axes2[0, 1].text(-0.1, 1.1, 'B', transform=axes2[0, 1].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

    # Lasso
    axes2[1, 0].scatter(y_test, results['Lasso']['y_pred'], alpha=0.5, color=palette[1])
    axes2[1, 0].plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'r--', lw=2)
    axes2[1, 0].set_title(f"Lasso Regression (MAE: {results['Lasso']['mae']:.2f}, R2: {results['Lasso']['r2']:.2f})", fontsize=14)
    axes2[1, 0].set_xlabel('Actual Age', fontsize=12)
    axes2[1, 0].set_ylabel('Predicted Age', fontsize=12)
    axes2[1, 0].text(-0.1, 1.1, 'C', transform=axes2[1, 0].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
    
    sns.residplot(x=y_test, y=results['Lasso']['y_pred'], ax=axes2[1, 1], color=palette[1])
    axes2[1, 1].set_title('Lasso Residuals', fontsize=14)
    axes2[1, 1].set_xlabel('Actual Age', fontsize=12)
    axes2[1, 1].set_ylabel('Residuals', fontsize=12)
    axes2[1, 1].text(-0.1, 1.1, 'D', transform=axes2[1, 1].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
    
    fig2.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig2.savefig('figure2_classical_models_refined.png', dpi=300)
    plt.close(fig2)

    # Figure 3: Gradient Boosting Models
    fig3, axes3 = plt.subplots(2, 2, figsize=(12, 10))
    fig3.suptitle('Figure 3: Gradient Boosting Models Performance', fontsize=20, fontweight='bold')
    
    # XGBoost
    axes3[0, 0].scatter(y_test, results['XGBoost']['y_pred'], alpha=0.5, color=palette[2])
    axes3[0, 0].plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'r--', lw=2)
    axes3[0, 0].set_title(f"XGBoost (MAE: {results['XGBoost']['mae']:.2f}, R2: {results['XGBoost']['r2']:.2f})", fontsize=14)
    axes3[0, 0].set_xlabel('Actual Age', fontsize=12)
    axes3[0, 0].set_ylabel('Predicted Age', fontsize=12)
    axes3[0, 0].text(-0.1, 1.1, 'A', transform=axes3[0, 0].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
    
    sns.residplot(x=y_test, y=results['XGBoost']['y_pred'], ax=axes3[0, 1], color=palette[2])
    axes3[0, 1].set_title('XGBoost Residuals', fontsize=14)
    axes3[0, 1].set_xlabel('Actual Age', fontsize=12)
    axes3[0, 1].set_ylabel('Residuals', fontsize=12)
    axes3[0, 1].text(-0.1, 1.1, 'B', transform=axes3[0, 1].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

    # LightGBM
    axes3[1, 0].scatter(y_test, results['LightGBM']['y_pred'], alpha=0.5, color=palette[3])
    axes3[1, 0].plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'r--', lw=2)
    axes3[1, 0].set_title(f"LightGBM (MAE: {results['LightGBM']['mae']:.2f}, R2: {results['LightGBM']['r2']:.2f})", fontsize=14)
    axes3[1, 0].set_xlabel('Actual Age', fontsize=12)
    axes3[1, 0].set_ylabel('Predicted Age', fontsize=12)
    axes3[1, 0].text(-0.1, 1.1, 'C', transform=axes3[1, 0].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
    
    sns.residplot(x=y_test, y=results['LightGBM']['y_pred'], ax=axes3[1, 1], color=palette[3])
    axes3[1, 1].set_title('LightGBM Residuals', fontsize=14)
    axes3[1, 1].set_xlabel('Actual Age', fontsize=12)
    axes3[1, 1].set_ylabel('Residuals', fontsize=12)
    axes3[1, 1].text(-0.1, 1.1, 'D', transform=axes3[1, 1].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
    
    fig3.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig3.savefig('figure3_boosting_models_refined.png', dpi=300)
    plt.close(fig3)

    # Figure 4: Model Comparison and Feature Importance
    fig4, axes4 = plt.subplots(1, 3, figsize=(18, 6))
    fig4.suptitle('Figure 4: Model Comparison and Feature Importance', fontsize=20, fontweight='bold')
    
    # Panel A: MAE Comparison
    mae_scores = {name: res['mae'] for name, res in results.items()}
    sns.barplot(x=list(mae_scores.keys()), y=list(mae_scores.values()), ax=axes4[0], palette="viridis")
    axes4[0].set_title('Mean Absolute Error (MAE)', fontsize=16)
    axes4[0].set_ylabel('MAE', fontsize=12)
    axes4[0].tick_params(axis='x', rotation=45)
    axes4[0].text(-0.1, 1.1, 'A', transform=axes4[0].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
    
    # Panel B: R-squared Comparison
    r2_scores = {name: res['r2'] for name, res in results.items()}
    sns.barplot(x=list(r2_scores.keys()), y=list(r2_scores.values()), ax=axes4[1], palette="viridis")
    axes4[1].set_title('R-squared', fontsize=16)
    axes4[1].set_ylabel('R-squared', fontsize=12)
    axes4[1].tick_params(axis='x', rotation=45)
    axes4[1].text(-0.1, 1.1, 'B', transform=axes4[1].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

    # Panel C: Feature Importance (for best model)
    best_model_name = max(r2_scores, key=r2_scores.get)
    best_model = results[best_model_name]['model']
    if hasattr(best_model, 'feature_importances_'):
        importances = best_model.feature_importances_
        feature_names = [f'PC{i+1}' for i in range(len(importances))]
        feature_importance_df = pd.DataFrame({'feature': feature_names, 'importance': importances})
        feature_importance_df = feature_importance_df.sort_values(by='importance', ascending=False).head(10)
        sns.barplot(x='importance', y='feature', data=feature_importance_df, ax=axes4[2], palette="viridis")
        axes4[2].set_title(f'Top 10 Feature Importances ({best_model_name})', fontsize=14)
        axes4[2].set_xlabel('Importance', fontsize=12)
        axes4[2].set_ylabel('Feature', fontsize=12)
    else:
        axes4[2].text(0.5, 0.5, 'Feature importance not available', ha='center', va='center')
    axes4[2].text(-0.1, 1.1, 'C', transform=axes4[2].transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

    fig4.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig4.savefig('figure4_comparison_refined.png', dpi=300)
    plt.close(fig4)
    print("Refined figures generated.")

if __name__ == '__main__':
    file_path = '/Users/vinayak/Documents/gemini-cli/data/archs4_human_counts_age_sex_organ.h5ad'
    X_pca, y, pca, df = load_and_preprocess_data(file_path)
    X_train, X_test, y_train, y_test = train_test_split(X_pca, y, test_size=0.2, random_state=42)
    results = train_and_evaluate_models(X_train, X_test, y_train, y_test)
    generate_refined_figures(results, pca, df, y_test)
    print("Analysis complete.")
