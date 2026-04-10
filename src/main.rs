use eframe::egui::{self, Color32, Pos2, Rect, Sense, Stroke, Vec2};
use crystal_sim::config::SimConfig;
use crystal_sim::kmc::Simulation;

const DEFAULT_CONFIG: &str = r##"{
    "particle_types": [
        {"radius": 0.5, "color": "#4488ff", "mu": 3.5},
        {"radius": 0.5, "color": "#ff6644", "mu": 3.5},
        {"radius": 0.5, "color": "#0cee44", "mu": 3.5}

    ],
    "epsilon": [[2.0, 6.0, 1.0], [6.0, 2.0, 1.0], [-10.0, -10.0, -10.0]],
    "delta": 0.05,
    "temperature": 6.0,
    "nu": 1.0,
    "seed": 42,
    "num_isolated_angles": 16,
    "relax_steps": 0,
    "relax_alpha": 0.01,
    "spring_k": 50.0
}"##;

// const DEFAULT_CONFIG: &str = r##"{
//     "particle_types": [
//         {"radius": 0.5, "color": "#4488ff", "mu": 3.5},
//         {"radius": 0.5, "color": "#ff6644", "mu": 3.5},
//     ],
//     "epsilon": [[-10.0, 1.0], [1.0, -10.0]],
//     "delta": 0.15,
//     "temperature": 1.0,
//     "nu": 1.0,
//     "seed": 42,
//     "num_isolated_angles": 16
// }"##;

fn parse_color(hex: &str) -> Color32 {
    let hex = hex.trim_start_matches('#');
    let r = u8::from_str_radix(&hex[0..2], 16).unwrap_or(128);
    let g = u8::from_str_radix(&hex[2..4], 16).unwrap_or(128);
    let b = u8::from_str_radix(&hex[4..6], 16).unwrap_or(128);
    Color32::from_rgb(r, g, b)
}

struct CrystalApp {
    sim: Simulation,
    config_string: String,
    running: bool,
    steps_per_frame: u32,
    zoom: f32,
    pan: Vec2,
    colors: Vec<Color32>,
    temperature: f32,
    mu: Vec<f32>,
    show_bonds: bool,
    show_candidates: bool,
    show_orientation: bool,
}

impl CrystalApp {
    fn new(config_string: String) -> Self {
        let mut config: SimConfig =
            serde_json::from_str(&config_string).expect("config invalid");
        config.init_cache();
        let colors = config.particle_types.iter().map(|t| parse_color(&t.color)).collect();
        let temperature = config.temperature;
        let mu = config.particle_types.iter().map(|t| t.mu).collect();
        let steps_per_frame = config.steps_per_frame;
        CrystalApp {
            sim: Simulation::new(config),
            config_string,
            running: true,
            steps_per_frame: steps_per_frame,
            zoom: 20.0,
            pan: Vec2::ZERO,
            colors,
            temperature,
            mu,
            show_bonds: true,
            show_candidates: false,
            show_orientation: false,
        }
    }

    fn sim_to_screen(&self, sim_x: f32, sim_y: f32, rect: Rect) -> Pos2 {
        let center = rect.center();
        Pos2::new(
            center.x + (sim_x + self.pan.x) * self.zoom,
            center.y - (sim_y + self.pan.y) * self.zoom,
        )
    }

    fn auto_center(&mut self) {
        if self.sim.particle_count() == 0 {
            return;
        }
        let mut cx = 0.0f32;
        let mut cy = 0.0f32;
        let mut n = 0;
        for p in self.sim.particle_grid.iter() {
            cx += p.pos.x;
            cy += p.pos.y;
            n += 1;
        }
        let n = n as f32;
        self.pan = Vec2::new(-(cx / n), -(cy / n));
    }
}

impl eframe::App for CrystalApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        if self.running {
            self.sim.step(self.steps_per_frame);
        }

        egui::SidePanel::left("controls").min_width(210.0).show(ctx, |ui| {
            ui.add_space(4.0);
            ui.heading("Crystal Simulator");
            ui.separator();

            ui.horizontal(|ui| {
                let label = if self.running { "Pause" } else { "Run" };
                if ui.button(label).clicked() {
                    self.running = !self.running;
                }
                if ui.button("1").clicked() {
                    self.sim.step(1);
                }
                if ui.button("10").clicked() {
                    self.sim.step(10);
                }
                if ui.button("100").clicked() {
                    self.sim.step(100);
                }
                if ui.button("Reset").clicked() {
                    let mut config: SimConfig =
                        serde_json::from_str(&self.config_string).expect("invalid config");
                    config.init_cache();
                    self.sim = Simulation::new(config);
                    self.pan = Vec2::ZERO;
                }
            });

            ui.separator();

            ui.label("Steps per frame");
            let mut spf = self.steps_per_frame as f64;
            if ui
                .add(egui::Slider::new(&mut spf, 1.0..=5000.0).logarithmic(true))
                .changed()
            {
                self.steps_per_frame = spf as u32;
            }

            ui.separator();

            // ui.label("Temperature");
            // if ui
            //     .add(egui::Slider::new(&mut self.temperature, 0.1..=5.0).step_by(0.05))
            //     .changed()
            // {
            //     self.sim.set_temperature(self.temperature);
            // }

            ui.separator();

            for (i, mu) in self.mu.iter_mut().enumerate() {
                let color = self.colors.get(i).copied().unwrap_or(Color32::WHITE);
                ui.horizontal(|ui| {
                    ui.colored_label(color, format!("mu type {}", i));
                });
                let mut m = *mu;
                // if ui.add(egui::Slider::new(&mut m, 0.0..=8.0).step_by(0.1)).changed() {
                //     *mu = m;
                //     self.sim.set_chemical_potential(i, m);
                // }
            }

            ui.separator();

            ui.checkbox(&mut self.show_bonds, "Show bonds");
            ui.checkbox(&mut self.show_candidates, "Show candidates");
            ui.checkbox(&mut self.show_orientation, "Show orientation");

            ui.separator();

            if ui.button("Auto-center").clicked() {
                self.auto_center();
            }

            ui.separator();

            ui.label(format!("Particles: {}", self.sim.particle_count()));
            ui.label(format!("Sim time:  {:.3}", self.sim.simulation_time()));
            ui.label(format!("Zoom:      {:.1}px/unit", self.zoom));

            ui.add_space(8.0);
            ui.small("Scroll: zoom");
            ui.small("Drag: pan");
        });

        egui::CentralPanel::default().show(ctx, |ui| {
            let (response, painter) =
                ui.allocate_painter(ui.available_size(), Sense::drag());
            let rect = response.rect;

            // Pan with drag
            if response.dragged() {
                let d = response.drag_delta();
                self.pan.x += d.x / self.zoom;
                self.pan.y -= d.y / self.zoom;
            }

            // Zoom with scroll wheel
            let scroll = ctx.input(|i| i.smooth_scroll_delta.y);
            if scroll != 0.0 {
                self.zoom = (self.zoom * (1.0 + scroll * 0.008)).clamp(1.0, 500.0);
            }

            // Dark background
            painter.rect_filled(rect, 0.0, Color32::from_rgb(12, 12, 22));

            // Bonds
            // if self.show_bonds {
            //     for i in 0..self.sim.particles.len() {
            //         let nbrs = self.sim.get_neighbors(i);
            //         for j in nbrs {
            //             if j <= i {
            //                 continue; // draw each bond once
            //             }
            //             let a = self.sim_to_screen(
            //                 self.sim.particles[i].pos.x,
            //                 self.sim.particles[i].pos.y,
            //                 rect,
            //             );
            //             let b = self.sim_to_screen(
            //                 self.sim.particles[j].pos.x,
            //                 self.sim.particles[j].pos.y,
            //                 rect,
            //             );
            //             painter.line_segment(
            //                 [a, b],
            //                 Stroke::new(1.2, Color32::from_rgba_unmultiplied(200, 200, 220, 70)),
            //             );
            //         }
            //     }
            // }

            // Attach-rate halos: ring around each particle scaled by its attach_rate_sum.
            if self.show_candidates {
                let max_rate = self.sim.particle_grid.cells.iter()
                    .flat_map(|c| c.attach_rates.iter().copied())
                    .fold(0.0f64, f64::max);
                let max_rate = if max_rate > 0.0 { max_rate } else { 1.0 };

                for cell in &self.sim.particle_grid.cells {
                    for (idx, p) in cell.particles.iter().enumerate() {
                        let rate = cell.attach_rates[idx];
                        if rate == 0.0 { continue; }
                        let pos = self.sim_to_screen(p.pos.x, p.pos.y, rect);
                        if pos.x < rect.left() - 50.0 || pos.x > rect.right() + 50.0
                            || pos.y < rect.top() - 50.0 || pos.y > rect.bottom() + 50.0
                        { continue; }
                        let t = (rate / max_rate) as f32;
                        let screen_r = (p.radius + 0.15) * self.zoom;
                        let alpha = ((0.2 + t * 0.6) * 255.0) as u8;
                        let color = self.colors.get(p.type_id).copied().unwrap_or(Color32::WHITE);
                        let [r, g, b, _] = color.to_array();
                        painter.circle_stroke(
                            pos,
                            screen_r,
                            Stroke::new(1.5 + t * 2.0, Color32::from_rgba_unmultiplied(r, g, b, alpha)),
                        );
                    }
                }
            }

            // Particles
            for particle in self.sim.particle_grid.iter() {
                let pos = self.sim_to_screen(particle.pos.x, particle.pos.y, rect);
                // Skip particles outside the visible rect with a margin
                if pos.x < rect.left() - 50.0
                    || pos.x > rect.right() + 50.0
                    || pos.y < rect.top() - 50.0
                    || pos.y > rect.bottom() + 50.0
                {
                    continue;
                }
                let screen_r = (particle.radius - self.sim.config.delta) * self.zoom;
                if screen_r < 0.5 {
                    continue;
                }
                let color = self
                    .colors
                    .get(particle.type_id)
                    .copied()
                    .unwrap_or(Color32::WHITE);
                painter.circle_filled(pos, screen_r, color);
                if screen_r > 3.0 {
                    painter.circle_stroke(
                        pos,
                        screen_r,
                        Stroke::new(0.8, Color32::from_rgba_unmultiplied(255, 255, 255, 50)),
                    );
                }

                // Orientation arrow: line from center toward orientation direction.
                if self.show_orientation && screen_r > 3.0 {
                    let ox = particle.orientation.x;
                    let oy = -particle.orientation.y; // flip for screen space
                    let tip = Pos2::new(pos.x + ox * screen_r * 0.75, pos.y + oy * screen_r * 0.75);
                    let arrow_color = Color32::from_rgba_unmultiplied(255, 255, 255, 200);
                    painter.line_segment([pos, tip], Stroke::new(1.5, arrow_color));
                    // Small arrowhead: two lines at ±30° from the tip
                    let head_len = screen_r * 0.25;
                    let angle = oy.atan2(ox);
                    for sign in [-1.0f32, 1.0] {
                        let a = angle + sign * std::f32::consts::PI * 5.0 / 6.0;
                        let h = Pos2::new(tip.x + a.cos() * head_len, tip.y + a.sin() * head_len);
                        painter.line_segment([tip, h], Stroke::new(1.5, arrow_color));
                    }
                }

                // Patch markers: draw small squares on the particle boundary.
                if let Some(type_def) = self.sim.config.particle_types.get(particle.type_id) {
                    let patch_half = (screen_r * 0.18).max(1.5);
                    for patch in &type_def.patches {
                        let n = crystal_sim::forces::patch_dir(particle.orientation, patch.position_cs);
                        let ux = n.x;
                        let uy = -n.y; // screen-space y is flipped
                        let vx = -uy;
                        let vy = ux;

                        let cx = pos.x + ux * screen_r;
                        let cy = pos.y + uy * screen_r;
                        let c = Pos2::new(cx, cy);

                        let p0 = Pos2::new(c.x - ux * patch_half - vx * patch_half, c.y - uy * patch_half - vy * patch_half);
                        let p1 = Pos2::new(c.x + ux * patch_half - vx * patch_half, c.y + uy * patch_half - vy * patch_half);
                        let p2 = Pos2::new(c.x + ux * patch_half + vx * patch_half, c.y + uy * patch_half + vy * patch_half);
                        let p3 = Pos2::new(c.x - ux * patch_half + vx * patch_half, c.y - uy * patch_half + vy * patch_half);

                        let patch_color = self.sim.config.patch_colors
                            .get(patch.patch_type_id)
                            .and_then(|c| c.as_deref())
                            .map(parse_color)
                            .unwrap_or(Color32::from_rgb(245, 245, 245));
                        painter.add(egui::Shape::convex_polygon(
                            vec![p0, p1, p2, p3],
                            patch_color,
                            Stroke::new(0.8, Color32::from_rgba_unmultiplied(20, 20, 20, 180)),
                        ));
                    }
                }
            }
        });

        // Keep animating
        ctx.request_repaint();
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let config_string = if args.len() > 1 {
        std::fs::read_to_string(&args[1]).expect("Failed to read config file")
    } else {
        DEFAULT_CONFIG.to_string()
    };

    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_title("Crystal Growth Simulator")
            .with_inner_size([1100.0, 750.0]),
        ..Default::default()
    };
    eframe::run_native(
        "Crystal Growth Simulator",
        options,
        Box::new(move |_cc| Ok(Box::new(CrystalApp::new(config_string)))),
    )
    .expect("failed to launch GUI");
}
