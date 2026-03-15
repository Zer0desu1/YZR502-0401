# 2-DOF Planar Robot Arm Simulation

Bu proje, 2 serbestlik dereceli (2-DOF) düzlemsel bir robot kolunun dinamik analizini ve ters dinamik kontrol simülasyonunu içermektedir. `robot_simulation.py` dosyası ile çalıştırılır.

## Özellikler

- **Dinamik Model**: Robotun atalet (mass), coriolis ve yerçekimi matrisleri hesaplanır.
- **Ters Dinamik Kontrol (Computed Torque Control)**: PID denetleyicisi ile robot kolunun istenen yörüngeyi takip etmesi sağlanır.
- **Yörünge Planlama**: Uç işlevci (end-effector) için görev uzayında (task space) dairesel bir yörünge oluşturulur.
- **Ters Kinematik (Inverse Kinematics)**: Görev uzayındaki konumlar, eklem açılarına dönüştürülür.
- **Simülasyon Çıktıları**: Simülasyon sonucunda aşağıdaki grafikler otomatik olarak proje klasörüne PNG formatında kaydedilir:
  - `fig_trajectory.png`: Uç işlevcinin istenen ve gerçekleşen dairesel yörüngesi.
  - `fig_joint_angles.png`: Her iki eklem için açı değişimleri.
  - `fig_errors.png`: Eklemlerin zaman bazlı yörünge takibi hata grafiği.
  - `fig_torques.png`: Her bir ekleme uygulanan kontrol torkları.
  - `fig_rms.png`: Eklemler için RMS takip hatası çubuk grafiği.
  - `fig_gain_analysis.png`: Pozisyon kazancı (Kp) ile RMS hata oranı ilişkisi grafiği.

## Kullanım

Eğer kendi ortamınızda test etmek istiyorsanız, gerekli Python kütüphanelerini yükleyin (`numpy`, `matplotlib`, `scipy` vb.) ve script dosyasını çalıştırın:

```bash
python robot_simulation.py
```
