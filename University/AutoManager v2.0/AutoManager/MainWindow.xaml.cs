using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using System.Windows.Threading;
using System.Drawing;
using System.Xml;
using System.Runtime.InteropServices;

using System.Windows.Forms;
using System.IO;
using System.Drawing.Imaging;

namespace AutoManager {
	public partial class MainWindow : Window {
		List<CppLib.ProcessInfo> Proc;
		DispatcherTimer timerStatusBar;
		List<ScriptItem> UserList;
		private delegate void UpdateProgressBarDelegate(System.Windows.DependencyProperty dp, Object value);
		CppLib.LowLevelKeyboardProcDelegate m_callback;
		Options options;
		bool Stop = false;

		public MainWindow() {
			InitializeComponent();
			options = LoadOptions();
			UpdatecbProgram();
			LoadcbScript();

			gvcCommand.Width = gvcProgram.Width = gvcTypeCommand.Width = window.Width / 3.3;

			sblScreen.Content = "X: " + SystemParameters.PrimaryScreenWidth + " Y: " + SystemParameters.PrimaryScreenHeight;
			System.Drawing.Point mousePosition = System.Windows.Forms.Cursor.Position;
			sblMousePoint.Content = "X: " + mousePosition.X + " Y: " + mousePosition.Y;
			timerStatusBar = new DispatcherTimer();
			timerStatusBar.Tick += new EventHandler(timerStatusBar_Tick);
			timerStatusBar.Interval = new TimeSpan(0, 0, 1);
			timerStatusBar.Start();
			UserList = new List<ScriptItem>();

			m_callback = LowLevelKeyboardHookProc;
			CppLib.SetHook(m_callback);
		}

		~MainWindow() {
			CppLib.Unhook();
		}

		private void timerStatusBar_Tick(object sender, EventArgs e) {
			System.Drawing.Point mousePosition = System.Windows.Forms.Cursor.Position;
			sblMousePoint.Content = "X: " + mousePosition.X + " Y: " + mousePosition.Y;
		}

		private void cbScript_SelectionChanged(object sender, SelectionChangedEventArgs e) {
			sblInfo.Content = "";
			if (cbScript.SelectedIndex != -1 && cbScript.SelectedIndex != cbScript.Items.Count-1) {
				LoadProg(cbScript.SelectedIndex);
				bSaveScript.Content = "Обновить сценарий";
				bDeleteScript.IsEnabled = true;
			}
			else {
				bSaveScript.Content = "Сохранить сценарий";
				bDeleteScript.IsEnabled = false;
				lvScript.Items.Clear();
				for (int i = 0; i < UserList.Count; i++)
					lvScript.Items.Add(UserList[i]);
			}
			bExeScript.IsEnabled = (lvScript.Items.Count != 0);
		}

		private void cbProgram_SelectionChanged(object sender, SelectionChangedEventArgs e) {
			bAddCommand.IsEnabled = false;
			if (cbProgram.SelectedItem != null) {
				if (cbProgram.SelectedItem.ToString() == "Системные команды") {
					lTypeCommand.Content = "Тип команды: ";
					cbTypeCommand.Items.Clear();
					cbTypeCommand.Items.Add("Кнопки мыши");
					cbTypeCommand.Items.Add("Перемещение мыши");
					cbTypeCommand.Items.Add("Ожидание");
					cbTypeCommand.Items.Add("Повтор");
				}
				else {
					lTypeCommand.Content = "Компонент: ";
					cbTypeCommand.Items.Clear();

					CppLib.ProcessInfo buff = Proc.Find(x => x.Name == cbProgram.SelectedItem.ToString());
					List<CppLib.ChildEll> ell = CppLib.GetNeedChildEll(buff.Handle, options.Ellements.ToArray());
					for (int i = 0; i < ell.Count; i++) {
						if (ell[i].Name == "")
							cbTypeCommand.Items.Add(i + "." + ell[i].ClassEll);
						else
							cbTypeCommand.Items.Add(i + "." +  ell[i].ClassEll + " \"" + ell[i].Name.Replace("&", "") + "\"");
					}

					cbTypeCommand.Items.Add(ell.Count + ".Контекстное меню");
				}
			}
		}

		private void cbTypeCommand_SelectionChanged(object sender, SelectionChangedEventArgs e) {
			cbCommand.Text = "";
			cbCommand.Items.Clear();
			bAddCommand.IsEnabled = cbCommand.IsEditable = false;
			if (cbTypeCommand.SelectedItem != null) {
				if (cbProgram.SelectedItem.ToString() == "Системные команды") {
					switch (cbTypeCommand.SelectedIndex) {
						case 0:
							cbCommand.IsEditable = false;
							cbCommand.Items.Add("Правая кнопка мыши");
							cbCommand.Items.Add("Левая кнопка мыши");
							cbCommand.Items.Add("Правая кнопка мыши вниз");
							cbCommand.Items.Add("Правая кнопка мыши вверх");
							cbCommand.Items.Add("Левая кнопка мыши вниз");
							cbCommand.Items.Add("Левая кнопка мыши вверх");
							break;
						case 1:
							cbCommand.IsEditable = true;
							break;
						case 2:
							cbCommand.IsEditable = true;
							for (Double i = 0.5; i <= 5; i = i + 0.5)
								cbCommand.Items.Add((i + " сек").Replace(',','.'));
							break;
						case 3:
							cbCommand.IsEditable = true;
							break;
					}
				}
				else {
					string buff = cbTypeCommand.SelectedItem.ToString().Split(' ')[0];
					buff = buff.Split('.')[1];
					switch (buff) {
						case "Edit":
							cbCommand.IsEditable = true;
							buff = "";
							if (cbTypeCommand.SelectedItem.ToString().Split(' ').Length > 1) {
								buff = cbTypeCommand.SelectedItem.ToString().Remove(0, cbTypeCommand.SelectedItem.ToString().IndexOf(' ') + 2);
								buff = buff.Remove(buff.Length - 1, 1);
							}
							cbCommand.Items.Add(buff);
							break;
						case "Button":
							cbCommand.IsEditable = false;
							cbCommand.Items.Add("Нажать кнопку");
							break;
						case "Контекстное":
							cbCommand.IsEditable = true;
							break;
					}
					cbCommand.SelectedIndex = 0;
				}
			}
		}

		private void cbCommand_SelectionChanged(object sender, SelectionChangedEventArgs e) {
			bAddCommand.IsEnabled = (cbProgram.SelectedIndex != -1 && cbTypeCommand.SelectedIndex != -1 && cbCommand.SelectedIndex != -1);
		}

		private void cbCommand_KeyUp(object sender, System.Windows.Input.KeyEventArgs e) {
			bAddCommand.IsEnabled = (cbProgram.SelectedIndex != -1 && cbTypeCommand.SelectedIndex != -1 && cbCommand.Text != "");
			if (e.Key == Key.Enter) {
				if (lvScript.SelectedIndex == -1 && bAddCommand.IsEnabled == true) {
					bAddCommand.RaiseEvent(new RoutedEventArgs(System.Windows.Controls.Button.ClickEvent));
				}
				else if(lvScript.SelectedIndex != -1 && bUpdateCommand.IsEnabled == true) {
					bUpdateCommand.RaiseEvent(new RoutedEventArgs(System.Windows.Controls.Button.ClickEvent));
				}
			}
		}

		private void lvScript_SelectionChanged(object sender, SelectionChangedEventArgs e) {
			bDeleteCommand.IsEnabled = bUpdateCommand.IsEnabled = (lvScript.SelectedIndex != -1);
			if(lvScript.SelectedIndex != -1) {
				ScriptItem item = (ScriptItem)lvScript.Items[lvScript.SelectedIndex];
				if(item.Program == "Системные команды") {
					cbProgram.SelectedIndex = cbProgram.Items.Count - 1;
					switch (item.TypeCommand) {
						case "Кнопки мыши":
							cbTypeCommand.SelectedIndex = 0;
							break;
						case "Перемещение мыши":
							cbTypeCommand.SelectedIndex = 1;
							break;
						case "Ожидание":
							cbTypeCommand.SelectedIndex = 2;
							break;
						case "Повтор":
							cbTypeCommand.SelectedIndex = 3;
							break;
					}
					cbCommand.Text = item.Command;
					bAddCommand.IsEnabled = true;
				}
			}
			else {
				cbProgram.SelectedIndex = cbProgram.Items.Count - 1;
				cbTypeCommand.SelectedIndex = cbCommand.SelectedIndex = -1;
			}
		}

		private void lvScript_MouseDown(object sender, MouseButtonEventArgs e) {
			lvScript.SelectedIndex = -1;
		}

		private void bUpdate_Click(object sender, RoutedEventArgs e) {
			UpdatecbProgram();
		}

		private void bAddCommand_Click(object sender, RoutedEventArgs e) {
			ScriptItem item = ValidScriptForm();

			if (item != null) {
				bExeScript.IsEnabled = true;

				if (lvScript.SelectedIndex == -1) {
					lvScript.Items.Add(item);
					UserList.Add(item);
				}
				else {
					lvScript.Items.Insert(lvScript.SelectedIndex + 1, item);
					lvScript.SelectedIndex = -1;
					if (cbScript.SelectedIndex == cbScript.Items.Count - 1) {
						UserList.Insert(lvScript.SelectedIndex + 1, item);
					}
				}
			}
		}

		private void bDeleteCommand_Click(object sender, RoutedEventArgs e) {
			if (cbScript.SelectedIndex == -1 || cbScript.SelectedIndex == cbScript.Items.Count - 1)
				UserList.RemoveAt(lvScript.SelectedIndex);
			lvScript.Items.RemoveAt(lvScript.SelectedIndex);
			bExeScript.IsEnabled = (lvScript.Items.Count != 0);
		}

		private void bSaveScript_Click(object sender, RoutedEventArgs e) {
			if (cbScript.SelectedIndex == -1 || cbScript.SelectedIndex == cbScript.Items.Count - 1)
				SaveProg();
			else
				UpdateScript(cbScript.SelectedIndex);
			LoadcbScript();
			UserList.Clear();
			lvScript.Items.Clear();
			cbTypeCommand.SelectedIndex = cbCommand.SelectedIndex = -1;
			bExeScript.IsEnabled = false;
		}

		private void bDeleteScript_Click(object sender, RoutedEventArgs e) {
			DeleteProg(cbScript.SelectedIndex);
			LoadcbScript();
			lvScript.Items.Clear();
		}
		private void bExeScript_Click(object sender, RoutedEventArgs e) {
			List<int[]> Repeat = new List<int[]>();
			int counter = options.CounterStart;
			int counter2 = 0;
			UpdatecbProgram();
			Stop = false;
			for (int i = 0; i < lvScript.Items.Count;i++) {
				lvScript.SelectedIndex = i;

				UpdateProgressBarDelegate updatePbDelegate = new UpdateProgressBarDelegate(lvScript.SetValue);
				Dispatcher.Invoke(updatePbDelegate, DispatcherPriority.Background, new object[] { System.Windows.Controls.ListView.SelectedIndexProperty, i });

				UpdateProgressBarDelegate updatePbDelegate2 = new UpdateProgressBarDelegate(sblInfo.SetValue);
				if(Repeat.Count > counter2)
					Dispatcher.Invoke(updatePbDelegate2, DispatcherPriority.Background, new object[] { System.Windows.Controls.Label.ContentProperty, "Осталось ещё " + Repeat[counter2][0] + " повторов" });

				ScriptItem item = new ScriptItem { Program = ((ScriptItem)lvScript.Items[i]).Program, TypeCommand = ((ScriptItem)lvScript.Items[i]).TypeCommand, Command = ((ScriptItem)lvScript.Items[i]).Command };

				if (Stop == true) {
					sblInfo.Content = "";
					Stop = false;
					lvScript.SelectedIndex = lvScript.Items.Count - 1;
					return;
				}

				while (true) {
					if (cbCounter.IsChecked == true && item.Command.Contains(options.CounterText) == true) {
						string buff2 = item.Command.Remove(item.Command.IndexOf(options.CounterText));
						buff2 = buff2 + counter + item.Command.Remove(0, item.Command.IndexOf(options.CounterText) + options.CounterText.Length);
						item.Command = buff2;
					}
					else
						break;
				}

				switch (item.TypeCommand) {
					case "Кнопки мыши":
						switch (item.Command) {
							case "Правая кнопка мыши":
								CppLib.MouseClickOnDesktop(MyMouseEvent.Mouse_Right_Click);
								break;
							case "Левая кнопка мыши":
								CppLib.MouseClickOnDesktop(MyMouseEvent.Mouse_Left_Click);
								break;
							case "Правая кнопка мыши вниз":
								CppLib.MouseClickOnDesktop(MyMouseEvent.Mouse_Right_Down);
								break;
							case "Правая кнопка мыши вверх":
								CppLib.MouseClickOnDesktop(MyMouseEvent.Mouse_Right_Up);
								break;
							case "Левая кнопка мыши вниз":
								CppLib.MouseClickOnDesktop(MyMouseEvent.Mouse_Left_Down);
								break;
							case "Левая кнопка мыши вверх":
								CppLib.MouseClickOnDesktop(MyMouseEvent.Mouse_Left_Up);
								break;
						}
						break;
					case "Перемещение мыши":
						System.Windows.Forms.Cursor.Position = new System.Drawing.Point(Convert.ToInt32(item.Command.Split(' ')[1]), Convert.ToInt32(item.Command.Split(' ')[3]));
						break;
					case "Ожидание":
						long tiks = DateTime.Now.Ticks + (long)(Convert.ToDouble(item.Command.Split(' ')[0].Replace('.',',')) * 10000000);
						while (true) {
							if (DateTime.Now.Ticks >= tiks)
								break;

							Dispatcher.Invoke(updatePbDelegate, DispatcherPriority.Background, new object[] { System.Windows.Controls.ListView.SelectedIndexProperty, i });

							if (Stop == true) {
								sblInfo.Content = "";
								Stop = false;
								lvScript.SelectedIndex = lvScript.Items.Count - 1;
								return;
							}
						}	
						break;
					case "Повтор":
						counter = counter + options.CounterStep;
						counter2 = 0;
						for (int j = 0; j < Repeat.Count + 1; j++) {
							if(j == Repeat.Count) {
								for(int k = 0; k < Repeat.Count;k++)
									Repeat[k][0] = Repeat[k][1];
								Repeat.Add(new int[3] { Convert.ToInt32(item.Command), Convert.ToInt32(item.Command), i });
								counter2 = j;
								break;
							}
							if (Repeat[j][2] == i) {
								counter2 = j;
								break;
							}
						}

						if(Repeat[counter2][0] != 0) {
							Repeat[counter2][0]--;
							i = -1;
							continue;
						}
						break;
					default:
						UpdatecbProgram();
						
						CppLib.ProcessInfo Process = Proc.Find(x => x.Name == item.Program);

						List<CppLib.ChildEll> ellements = CppLib.GetNeedChildEll(Process.Handle, options.Ellements.ToArray());
						string spliteType = item.TypeCommand.Split(' ')[0];
						string[] splitCommand = item.Command.Split(' ');
						string number = spliteType.Remove(spliteType.IndexOf('.'),spliteType.Length-1);
						string name = spliteType.Remove(0, spliteType.IndexOf('.')+1);
						if (name == "Контекстное") {
							Bitmap screenAfterClick, screenBeforClick;
							Graphics graphScreenAfterClick, graphScreenBeforClick;
							System.Drawing.Color[] colorsScreenAfterClick, colorsScreenBeforClick;

							screenAfterClick = new Bitmap(Process.Border.right - Process.Border.left, Process.Border.bottom - Process.Border.top);
							graphScreenAfterClick = Graphics.FromImage(screenAfterClick);

							screenBeforClick = new Bitmap(Process.Border.right - Process.Border.left, Process.Border.bottom - Process.Border.top);
							graphScreenBeforClick = Graphics.FromImage(screenBeforClick);


							graphScreenBeforClick.CopyFromScreen(new System.Drawing.Point(Process.Border.left, Process.Border.top), System.Drawing.Point.Empty, new System.Drawing.Size(Process.Border.right - Process.Border.left, Process.Border.bottom - Process.Border.top));

							System.Windows.Forms.Cursor.Position = new System.Drawing.Point(Convert.ToInt32(splitCommand[1]), Convert.ToInt32(splitCommand[3]));
							CppLib.MouseClickOnDesktop(MyMouseEvent.Mouse_Right_Click);

							long tiks2 = DateTime.Now.Ticks + (long)(Convert.ToDouble(0.5 * 10000000));

							colorsScreenBeforClick = GetColorArr(screenBeforClick);

							while (true) {
								if (DateTime.Now.Ticks >= tiks2)
									break;

								if (Stop == true) {
									sblInfo.Content = "";
									Stop = false;
									lvScript.SelectedIndex = lvScript.Items.Count - 1;
									return;
								}
							}

							graphScreenAfterClick.CopyFromScreen(new System.Drawing.Point(Process.Border.left, Process.Border.top), System.Drawing.Point.Empty, new System.Drawing.Size(Process.Border.right - Process.Border.left, Process.Border.bottom - Process.Border.top));
							colorsScreenAfterClick = GetColorArr(screenAfterClick);

							for (int y = 0; y < screenBeforClick.Height; y++) {
								for (int x = 0; x < screenBeforClick.Width; x++) {
									screenBeforClick.SetPixel(x, y, colorsScreenBeforClick[x + y * screenBeforClick.Width]);
								}
							}
							image1.Source = BitmapToImageSource(screenBeforClick);

							for (int y = 0; y < screenAfterClick.Height; y++) {
								for (int x = 0; x < screenAfterClick.Width; x++) {
									screenAfterClick.SetPixel(x, y, colorsScreenAfterClick[x + y * screenBeforClick.Width]);
								}
							}
							image.Source = BitmapToImageSource(screenAfterClick);
						}
						else {
							for (int j = 0; j < ellements.Count; j++) {
								if (name == "Edit" && j == Convert.ToInt32(number)) {
									for (int k = 0; k < cbTypeCommand.Items.Count; k++) {
										if (number + ".Edit" == cbTypeCommand.Items[k].ToString().Split(' ')[0]) {
											cbTypeCommand.Items[k] = number + ".Edit \"" + item.Command + "\"";
											break;
										}
									}
									CppLib.SetTitleText(ellements[j].Handle, item.Command);
									break;
								}
								if (name == "Button" && j == Convert.ToInt32(number)) {
									CppLib.ClickButton(ellements[j].Handle);
									break;
								}
							}
						}	
						break;
				}
			}
			sblInfo.Content = "";
		}

		private void bUpdateCommand_Click(object sender, RoutedEventArgs e) {
			ScriptItem item = ValidScriptForm();
			if (item != null) {
				lvScript.Items[lvScript.SelectedIndex] = item;
				lvScript.SelectedIndex = -1;
			}
		}

		private void bSaveOption_Click(object sender, RoutedEventArgs e) {
			SaveOption();
			cbEllements.SelectedIndex = cbChildWindow.SelectedIndex = -1;
			bAddChildWindow.IsEnabled = (cbChildWindow.Text.Length != 0);
			bAddElements.IsEnabled = (cbEllements.Text.Length != 0);
		}

		private void cbCounter_Checked(object sender, RoutedEventArgs e) {
			if (tbCounter != null) {
				tbCounter.IsEnabled = !tbCounter.IsEnabled;
				tbCounterStep.IsEnabled = !tbCounterStep.IsEnabled;
				tbCounterStart.IsEnabled = !tbCounterStart.IsEnabled;
				ValidCounterForm();
			}
		}

		private void cbChildWindow_SelectionChanged(object sender, SelectionChangedEventArgs e) {
			if (cbChildWindow.SelectedIndex == -1)
				bAddChildWindow.Content = "Добавить запись";
			else {
				bAddChildWindow.Content = "Удалить запись";
				bAddChildWindow.IsEnabled = true;
			}
		}

		private void cbEllements_SelectionChanged(object sender, SelectionChangedEventArgs e) {
			if (cbEllements.SelectedIndex == -1)
				bAddElements.Content = "Добавить запись";
			else {
				bAddElements.Content = "Удалить запись";
				bAddElements.IsEnabled = true;
			}
		}

		private void bAddChildWindow_Click(object sender, RoutedEventArgs e) {
			bAddChildWindow.IsEnabled = false;
			if (bAddChildWindow.Content.ToString() == "Добавить запись") {
				cbChildWindow.Items.Add(cbChildWindow.Text);
				options.ChildWindows.Add(cbChildWindow.Text);
				cbChildWindow.Text = "";
			}
			else {
				for (int i = 0; i < options.ChildWindows.Count; i++)
					if (options.ChildWindows[i].ToString() == cbChildWindow.Items[cbChildWindow.SelectedIndex].ToString()) {
						options.ChildWindows.RemoveAt(i);
						break;
					}
				cbChildWindow.Items.RemoveAt(cbChildWindow.SelectedIndex);
				cbChildWindow.SelectedIndex = -1;
			}
		}

		private void bAddElements_Click(object sender, RoutedEventArgs e) {
			bAddElements.IsEnabled = false;
			if (bAddElements.Content.ToString() == "Добавить запись") {
				cbEllements.Items.Add(cbEllements.Text);
				options.Ellements.Add(cbEllements.Text);
				cbEllements.Text = "";
			}
			else {
				for (int i = 0; i < options.Ellements.Count; i++)
					if (options.Ellements[i].ToString() == cbEllements.Items[cbEllements.SelectedIndex].ToString()) {
						options.Ellements.RemoveAt(i);
						break;
					}
				cbEllements.Items.RemoveAt(cbEllements.SelectedIndex);
				cbEllements.SelectedIndex = -1;
			}
		}

		private void cbChildWindow_KeyUp(object sender, System.Windows.Input.KeyEventArgs e) {
			bAddChildWindow.IsEnabled = !(cbChildWindow.SelectedIndex == -1 && cbChildWindow.Text.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries).Length == 0);
			if(e.Key == Key.Enter && bAddChildWindow.IsEnabled == true) {
				bAddChildWindow.RaiseEvent(new RoutedEventArgs(System.Windows.Controls.Button.ClickEvent));
			}
		}

		private void cbEllements_KeyUp(object sender, System.Windows.Input.KeyEventArgs e) {
			bAddElements.IsEnabled = !(cbEllements.SelectedIndex == -1 && cbEllements.Text.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries).Length == 0);
			if (e.Key == Key.Enter && bAddElements.IsEnabled == true) {
				bAddElements.RaiseEvent(new RoutedEventArgs(System.Windows.Controls.Button.ClickEvent));
			}
		}

		private void tbCounter_TextChanged(object sender, TextChangedEventArgs e) {
			ValidCounterForm();
		}

		private void tbCounterStart_TextChanged(object sender, TextChangedEventArgs e) {
			ValidCounterForm();
		}

		private void tbCounterStep_TextChanged(object sender, TextChangedEventArgs e) {
			ValidCounterForm();
		}

		private void button_Click(object sender, RoutedEventArgs e) {
			CppLib.Rect prob = new CppLib.Rect { bottom = 200, left = 100, right = 200, top = 100 };
			System.Drawing.Point prob2 = new System.Drawing.Point(450, 500); 

			Bitmap screenAfterClick, screenBeforClick;
			Graphics graphScreenAfterClick, graphScreenBeforClick;
			System.Drawing.Color[] colorsScreenAfterClick, colorsScreenBeforClick;

			screenAfterClick = new Bitmap(prob.right - prob.left, prob.bottom - prob.top);
			graphScreenAfterClick = Graphics.FromImage(screenAfterClick);

			screenBeforClick = new Bitmap(prob.right - prob.left, prob.bottom - prob.top);
			graphScreenBeforClick = Graphics.FromImage(screenBeforClick);


			graphScreenBeforClick.CopyFromScreen(new System.Drawing.Point(prob.left,prob.top), System.Drawing.Point.Empty, new System.Drawing.Size(prob.right-prob.left,prob.bottom-prob.top));

			System.Windows.Forms.Cursor.Position = new System.Drawing.Point(prob2.X, prob2.Y);
			CppLib.MouseClickOnDesktop(MyMouseEvent.Mouse_Right_Click);

			long tiks = DateTime.Now.Ticks + (long)(Convert.ToDouble(0.5 * 10000000));

			colorsScreenBeforClick = GetColorArr(screenBeforClick);

			while (true) {
				if (DateTime.Now.Ticks >= tiks)
					break;

				if (Stop == true) {
					sblInfo.Content = "";
					Stop = false;
					lvScript.SelectedIndex = lvScript.Items.Count - 1;
					return;
				}
			}

			graphScreenAfterClick.CopyFromScreen(new System.Drawing.Point(prob.left, prob.top), System.Drawing.Point.Empty, new System.Drawing.Size(prob.right - prob.left, prob.bottom - prob.top));
			colorsScreenAfterClick = GetColorArr(screenAfterClick);

			for (int y = 0; y < screenBeforClick.Height; y++) {
				for (int x = 0; x < screenBeforClick.Width; x++) {
					screenBeforClick.SetPixel(x, y, colorsScreenBeforClick[x + y * 100]);
				}
			}
			image1.Source = BitmapToImageSource(screenBeforClick);

			for (int y = 0; y < screenAfterClick.Height; y++) {
				for (int x = 0; x < screenAfterClick.Width; x++) {
					screenAfterClick.SetPixel(x, y, colorsScreenAfterClick[x + y * 100]);
				}
			}
			image.Source = BitmapToImageSource(screenAfterClick);
		}
	}
}