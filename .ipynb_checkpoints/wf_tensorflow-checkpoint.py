import uproot
import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc

# ----------------------
# 1. Load Data from ROOT
# ----------------------
print("Loading data from ROOT file...")
file = uproot.open("https://zenodo.org/records/16903129/files/wf_array.root?download=1")

tree = file["wf_array"]

wf_n_gamma = tree["wf_n_gamma"].array(library="np")  # shape: (N, 512)
wf_gamma = tree["wf_gamma"].array(library="np")      # shape: (N, 512)

print(f"Neutron waveforms: {wf_n_gamma.shape}, Gamma waveforms: {wf_gamma.shape}")

# Combine and label
X = np.vstack([wf_n_gamma, wf_gamma])  # (2N, 512)
y = np.hstack([np.ones(len(wf_n_gamma)), np.zeros(len(wf_gamma))])  # 1=neutron, 0=gamma

# Normalize and reshape for Conv1D
X = X / np.max(X)
X = X[..., np.newaxis]  # shape: (samples, 512, 1)

# ----------------------
# 2. Train/Test Split
# ----------------------
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

print(f"Train size: {X_train.shape}, Test size: {X_test.shape}")

# ----------------------
# 3. Build CNN Model
# ----------------------
model = tf.keras.Sequential([
    tf.keras.Input(shape=(512, 1)),
    tf.keras.layers.Conv1D(32, 5, activation='relu'),
    tf.keras.layers.MaxPooling1D(2),
    tf.keras.layers.Conv1D(64, 5, activation='relu'),
    tf.keras.layers.MaxPooling1D(2),
    tf.keras.layers.Flatten(),
    tf.keras.layers.Dense(128, activation='relu'),
    tf.keras.layers.Dropout(0.3),
    tf.keras.layers.Dense(1, activation='sigmoid')
])

model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])
model.summary()

# ----------------------
# 4. Train Model
# ----------------------
history = model.fit(X_train, y_train,
                    validation_split=0.2,
                    epochs=20,
                    batch_size=128,
                    verbose=1)

# ----------------------
# 5. Evaluate
# ----------------------
test_loss, test_acc = model.evaluate(X_test, y_test)
print(f"Test Accuracy: {test_acc*100:.2f}%")

# ----------------------
# 6. Visualizations
# ----------------------
# Accuracy and Loss curves
plt.figure(figsize=(10,4))
plt.subplot(1,2,1)
plt.plot(history.history['accuracy'], label='Train Acc')
plt.plot(history.history['val_accuracy'], label='Val Acc')
plt.xlabel('Epoch'); plt.ylabel('Accuracy'); plt.legend(); plt.title('Accuracy')

plt.subplot(1,2,2)
plt.plot(history.history['loss'], label='Train Loss')
plt.plot(history.history['val_loss'], label='Val Loss')
plt.xlabel('Epoch'); plt.ylabel('Loss'); plt.legend(); plt.title('Loss')
plt.tight_layout()
#plt.show()
plt.savefig("./output/cnn_training_curves.png")

# ROC Curve
y_pred_prob = model.predict(X_test).ravel()
fpr, tpr, _ = roc_curve(y_test, y_pred_prob)
roc_auc = auc(fpr, tpr)

plt.figure()
plt.plot(fpr, tpr, label=f"ROC Curve (AUC = {roc_auc:.2f})")
plt.plot([0,1], [0,1], 'k--')
plt.xlabel('False Positive Rate'); plt.ylabel('True Positive Rate')
plt.title('ROC Curve'); plt.legend()
#plt.show()
plt.savefig("./output/cnn_roc_curve.png")

# ----------------------
# 7. Save Model
# ----------------------
model.save("./output/cnn_waveform_classifier.h5")
print("Model saved as cnn_waveform_classifier.h5 in ./output directory")
print("Open the .png graphs in ./output with eog")