using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Media.Imaging;
using System.Diagnostics;
using System.Windows.Automation;
using System.Windows.Controls;
using System.Windows.Controls.Primitives;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Ink;
using System.Windows.Input;
using System.Windows.Markup;
using System.Windows.Media;
using System.Windows.Media.Animation;
using System.Windows.Media.Effects;
using System.Windows.Media.Media3D;
using System.Windows.Media.TextFormatting;
using System.Windows.Navigation;
using System.Windows.Shapes;
using System.Windows.Shell;
using System.Globalization;
using System.Windows.Forms;
using System.IO;
using System.Threading;
using System.Windows.Threading;

namespace prob3 {
	class PlayerForm {
		Image Form;
		MainWindow MainForm;
		DrawingGroup Player;
		int screenWidth = 500, screenHeight = 425, Index = 0;
		MyLabel Name, PathDerectory, Artists, Number, Album;
		MyButton Play, Next, Previous, Open;
		MyForm PlayerMyForm;
		MyProgressBar Bar;
		MyImage img;
		Rect[] ProgressBorder = { new Rect(10, 340, 484, 25), new Rect(0, 0, 28, 17) };
		Rect[] FormBorders = { new Rect(0, 0, 500, 30), new Rect(0, 40, 0, 0), new Rect(0, 5, 20, 20), new Rect(0, 5, 20, 20), new Rect(0, 5, 20, 20) };
		MediaPlayer mediaPlayer = new MediaPlayer();
		FolderBrowserDialog OpenDialog = new FolderBrowserDialog();
		List<string> MP3;
		TagLib.File FileInfo;
		BitmapImage bitmap;
		Thread thread;
		double animBar = 0.01;

		public PlayerForm(Image Form, MainWindow MainForm) {
			this.MainForm = MainForm;
			this.Form = Form;
			Player = new DrawingGroup();
			MP3 = new List<string>();
			bitmap = null;
			thread = new Thread(AnimationBar);
			thread.Start();
			//OpenDialog.Filter = "Audio files(*.mp3;*.wav) | *.mp3;*.wav";

			PlayerMyForm = new MyForm("Player", FormBorders, new Point(), screenWidth, screenHeight, Form, MainForm, Player);
			PlayerMyForm.EventFormClose = Form_Close;

			Name = new MyLabel("Name", " ", 14, 170, 85, Player);
			Name.WidthBody = 205;
			Name.WidthHeader = 95;
			Artists = new MyLabel("Artists", " ", 14, 170, 85 + Name.Height + 12, Player);
			Artists.WidthBody = 205;
			Artists.WidthHeader = 95;
			Album = new MyLabel("Album", " ", 14, 170, 85 + Name.Height * 2 + 24, Player);
			Album.WidthBody = 205;
			Album.WidthHeader = 95;
			Number = new MyLabel("Track Number", "1/1", 14, 170, 85 + Name.Height * 3 + 36, Player);
			Number.WidthBody = 205;
			Number.WidthHeader = 95;
			PathDerectory = new MyLabel("Path", " ", 14, 10, 245, Player);
			PathDerectory.WidthBody = 310;

			Previous = new MyButton(new Rect(30, 375, 115, 35), "Previous", 14, Player);
			Previous.Enabled = false;
			Previous.EventButtonKeyUp = PreviousKeyUp;
			Play = new MyButton(new Rect(195, 375, 115, 35), "Play", 14, Player);
			Play.Enabled = false;
			Play.EventButtonKeyUp = PlayKeyUp;
			Next = new MyButton(new Rect(360, 375, 115, 35), "Next", 14, Player);
			Next.Enabled = false;
			Next.EventButtonKeyUp = NextKeyUp;
			Open = new MyButton(new Rect(370, 240, 115, 35), "Open", 14, Player);
			Open.EventButtonKeyUp = OpenKeyUp;

			Bar = new MyProgressBar(ProgressBorder, Player);

			img = new MyImage(bitmap, new Rect(10, 85, 150, 150), Player);

			Form.Height = screenHeight;
			Form.Width = screenWidth;
			Form.Source = new DrawingImage(Player);

			mediaPlayer.MediaEnded += new EventHandler(player_MediaEnded);
		}

		public void Player_MouseDownKey(Point point) {
			PlayerMyForm.Form_MouseDownKey(point);
			Play.Button_MouseKeyDown(point);
			Next.Button_MouseKeyDown(point);
			Previous.Button_MouseKeyDown(point);
			Open.Button_MouseKeyDown(point);
		}

		public void Player_MouseUpKey(Point point) {
			PlayerMyForm.Form__MouseUpKey(point);
			Play.Button_MouseKeyUp(point);
			Next.Button_MouseKeyUp(point);
			Previous.Button_MouseKeyUp(point);
			Open.Button_MouseKeyUp(point);
		}

		public void Player_MouseMoveKey(Point point) {
			PlayerMyForm.Form__MouseMoveKey(point);
			Play.Button_MouseKeyMove(point);
			Next.Button_MouseKeyMove(point);
			Previous.Button_MouseKeyMove(point);
			Open.Button_MouseKeyMove(point);
		}

		void OpenKeyUp() {
			if (OpenDialog.ShowDialog() == DialogResult.OK) {
				PathDerectory.Text = OpenDialog.SelectedPath;
				MP3.Clear();
				string[] Files = Directory.GetFiles(OpenDialog.SelectedPath);
				for (int i = 0; i < Files.Length; i++) {
					string File = Files[i];
					if (File.Substring(File.Length - 3).ToLower() == "mp3") {
						MP3.Add(Files[i]);
					}
				}
				if (MP3.Count >= 1) {
					Play.Enabled = true;
					Index = 0;
					InitSong();
					if (MP3.Count >= 1) {
						Previous.Enabled = true;
						Next.Enabled = true;
					}
				}
			}
		}

		void PlayKeyUp() {
			if (Play.Text == "Play") {
				Play.Text = "Pause";
				mediaPlayer.Play();
			}
			else {
				Play.Text = "Play";
				mediaPlayer.Pause();
			}
		}

		void NextKeyUp() {
			if (Index == MP3.Count - 1) {
				Index = 0;
			}
			else {
				Index++;
			}
			InitSong();
		}

		void PreviousKeyUp() {
			if (Index == 0) {
				Index = MP3.Count - 1;
			}
			else {
				Index--;
			}
			InitSong();
		}

		void Form_Close() {
			thread.Abort();
		}

		private void player_MediaEnded(Object sender, EventArgs e) {
			NextKeyUp();
		}

		void InitSong() {
			FileInfo = TagLib.File.Create(MP3[Index]);
			Bar.Max = FileInfo.Properties.Duration.TotalSeconds;
			Name.Text = (FileInfo.Tag.Title != null) ? FileInfo.Tag.Title : " ";
			Artists.Text = (FileInfo.Tag.FirstAlbumArtist != null) ? FileInfo.Tag.FirstAlbumArtist : " ";
			Album.Text = (FileInfo.Tag.Album != null) ? FileInfo.Tag.Album : " ";
			Number.Text = (Index + 1).ToString() + "/" + MP3.Count.ToString();
			mediaPlayer.Open(new Uri(MP3[Index], UriKind.Absolute));
			Play.Text = "Pause";
			mediaPlayer.Play();
			bitmap = null;
			if (FileInfo.Tag.Pictures.Length != 0) {
				TagLib.IPicture pic = FileInfo.Tag.Pictures[0];
				MemoryStream ms = new MemoryStream(pic.Data.Data);
				ms.Seek(0, SeekOrigin.Begin);

				bitmap = new BitmapImage();
				bitmap.BeginInit();
				bitmap.StreamSource = ms;
				bitmap.EndInit();
			}
			img.Picture = bitmap;
		}

		private void AnimationBar() {
			while (true) {
				Thread.Sleep(TimeSpan.FromMilliseconds(80));
				MainForm.Dispatcher.BeginInvoke(DispatcherPriority.Normal,
					(ThreadStart)delegate () {
						if (Play.Text == "Pause") {
							if (Bar.BarOpacity >= 1)
								animBar = -0.05;
							if (Bar.BarOpacity <= 0.15)
								animBar = 0.05;
							Bar.BarOpacity = Bar.BarOpacity + animBar;
							Bar.Text = mediaPlayer.Position.TotalSeconds.ToString() + " / " + Bar.Max;
							Bar.Progress = mediaPlayer.Position.TotalSeconds;
						}
						if(Play.Text == "Play") {
							if(Bar.BarOpacity != 1) {
								Bar.BarOpacity = 1;
								Bar.Progress = mediaPlayer.Position.TotalSeconds;
							}
						}
					}
			);
			}
		}
	}
}
